"""Real-time three-array SRP-PHAT localization.

Example:
  python realtime_localizer.py --list-devices
  python realtime_localizer.py --devices 1 2 3 --channels 8 16 16 \
      --positions microphone_positions.csv

Keep three_array_localizer.py in the same folder. During calibration, play white
noise from a speaker equidistant from the selected sync microphone on each array.
Stop with Ctrl+C.
"""

from __future__ import annotations

import argparse
import queue
import time

import numpy as np
import sounddevice as sd

from three_array_localizer import (
    SyncResult,
    align_array,
    estimate_sync,
    SRPPhatLocalizer,
)


class LiveBuffer:
    """Sliding window over the most recent `capacity` samples.

    Storage is twice the capacity so `data` is always one contiguous view;
    the window is shifted back only once per `capacity` samples.
    """

    def __init__(self, channels: int, capacity: int) -> None:
        self._storage = np.zeros((2 * capacity, channels), dtype=np.float32)
        self._start = 0
        self._stop = 0
        self.base = 0
        self.capacity = capacity

    @property
    def data(self) -> np.ndarray:
        return self._storage[self._start : self._stop]

    @property
    def end(self) -> int:
        return self.base + self._stop - self._start

    def append(self, block: np.ndarray) -> None:
        n = len(block)
        if self._stop + n > len(self._storage):
            keep = self._storage[self._stop - self.capacity + n : self._stop].copy()
            self.base += self._stop - self._start - len(keep)
            self._storage[: len(keep)] = keep
            self._start, self._stop = 0, len(keep)
        self._storage[self._stop : self._stop + n] = block
        self._stop += n
        excess = self._stop - self._start - self.capacity
        if excess > 0:
            self._start += excess
            self.base += excess


class LivePlot:
    """Top (x-y) and side (y-z) views of the microphones and source track."""

    def __init__(self, positions: np.ndarray, channels, bounds: np.ndarray) -> None:
        import matplotlib.pyplot as plt

        self.plt = plt
        plt.ion()
        self.fig, (self.top, self.side) = plt.subplots(1, 2, figsize=(11, 5))
        names = ["UMA-8", "UMA-16 L", "UMA-16 R"]
        edges = np.cumsum([0, *channels])
        for name, lo, hi in zip(names, edges[:-1], edges[1:]):
            p = positions[lo:hi]
            self.top.scatter(p[:, 0], p[:, 1], s=8, label=name)
            self.side.scatter(p[:, 1], p[:, 2], s=8, label=name)
        for ax, (i, j), labels in (
            (self.top, (0, 1), ("x [m] (right)", "y [m] (front)")),
            (self.side, (1, 2), ("y [m] (front)", "z [m] (up)")),
        ):
            ax.set_xlim(*bounds[i])
            ax.set_ylim(*bounds[j])
            ax.set_xlabel(labels[0])
            ax.set_ylabel(labels[1])
            ax.set_aspect("equal", adjustable="box")
            ax.grid(alpha=0.3)
        lo = np.minimum(bounds[:, 0], positions.min(axis=0) - 0.1)
        hi = np.maximum(bounds[:, 1], positions.max(axis=0) + 0.1)
        self.top.set_xlim(lo[0], hi[0])
        self.top.set_ylim(lo[1], hi[1])
        self.side.set_xlim(lo[1], hi[1])
        self.side.set_ylim(lo[2], hi[2])
        self.top.set_title("Top view")
        self.side.set_title("Side view")
        self.top.legend(loc="lower right", fontsize=8)
        (self.top_trail,) = self.top.plot([], [], "-", color="0.6", lw=1)
        (self.side_trail,) = self.side.plot([], [], "-", color="0.6", lw=1)
        (self.top_dot,) = self.top.plot([], [], "o", color="red", ms=12)
        (self.side_dot,) = self.side.plot([], [], "o", color="red", ms=12)
        self.track: list[np.ndarray] = []
        self.fig.tight_layout()
        self.pause(0.05)

    def update(self, position: np.ndarray, score: float, level: float,
               clear: bool) -> None:
        """Move the dot every frame: red when clear, grey when the score is low."""
        if clear:
            self.track = (self.track + [position])[-20:]
            t = np.asarray(self.track)
            self.top_trail.set_data(t[:, 0], t[:, 1])
            self.side_trail.set_data(t[:, 1], t[:, 2])
        x, y, z = position
        for dot, xy in ((self.top_dot, ([x], [y])), (self.side_dot, ([y], [z]))):
            dot.set_data(*xy)
            dot.set_color("red" if clear else "0.6")
            dot.set_markersize(12 if clear else 8)
        status = "" if clear else "  (weak - no clear source)"
        self.fig.suptitle(
            f"Source x={x:+.2f}  y={y:+.2f}  z={z:+.2f} m   "
            f"score={score:.2f}  level={level:.0f} dBFS{status}"
        )

    def message(self, text: str) -> None:
        self.fig.suptitle(text)

    def pause(self, seconds: float) -> None:
        self.plt.pause(seconds)


def drain_queues(queues, buffers, sync_record=None, sync_channels=None) -> None:
    """Move captured blocks into the buffers.

    During calibration, sync_record collects the full sync channel of each array,
    since the live buffers only keep the most recent seconds.
    """
    for i, (q, buffer) in enumerate(zip(queues, buffers)):
        while True:
            try:
                block = q.get_nowait()
            except queue.Empty:
                break
            buffer.append(block)
            if sync_record is not None:
                sync_record[i].append(block[:, sync_channels[i]].copy())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list-devices", action="store_true")
    parser.add_argument("--devices", type=int, nargs=3, metavar=("A0", "A1", "A2"))
    parser.add_argument("--channels", type=int, nargs=3, metavar=("N0", "N1", "N2"))
    parser.add_argument("--positions", help="CSV: one x,y,z row per captured channel")
    parser.add_argument("--sync-channels", type=int, nargs=3, default=[0, 0, 0])
    parser.add_argument("--fs", type=int, default=48000)
    parser.add_argument("--calibration-seconds", type=float, default=45.0)
    parser.add_argument("--frame-ms", type=int, default=200)
    parser.add_argument("--update-ms", type=int, default=500)
    parser.add_argument("--coarse-step", type=float, default=0.10)
    parser.add_argument("--refinements", type=int, default=3)
    parser.add_argument("--min-score", type=float, default=0.25,
                        help="estimates below this SRP score are treated as no clear source")
    parser.add_argument("--plot", action="store_true",
                        help="show a live top/side view of the source position")
    parser.add_argument(
        "--bounds", type=float, nargs=6, default=[-1, 1, -1, 1, 0.05, 2],
        metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"),
    )
    args = parser.parse_args()

    if args.list_devices:
        print(sd.query_devices())
        return
    if not args.devices or not args.channels or not args.positions:
        parser.error("use --devices, --channels and --positions, or --list-devices")
    if args.calibration_seconds < 5:
        parser.error("use at least 5 calibration seconds")

    # Skip an optional "x,y,z" header row.
    with open(args.positions) as f:
        first = f.readline().split(",")[0].strip()
    try:
        float(first)
        header_rows = 0
    except ValueError:
        header_rows = 1
    positions = np.loadtxt(
        args.positions, delimiter=",", ndmin=2, skiprows=header_rows
    )
    total_channels = sum(args.channels)
    if positions.shape != (total_channels, 3):
        parser.error(
            f"{args.positions} must contain {total_channels} x,y,z rows (one per "
            f"channel), found {positions.shape[0] if positions.size else 0}. "
            "Run make_positions.py to generate it."
        )
    for device, channels in zip(args.devices, args.channels):
        info = sd.query_devices(device)
        if info["max_input_channels"] < channels:
            raise ValueError(f"Device {device} does not provide {channels} inputs")
    for channel, channels in zip(args.sync_channels, args.channels):
        if not 0 <= channel < channels:
            raise ValueError("A sync channel is outside its array's channel range")

    fs = args.fs
    frame_samples = int(args.frame_ms * fs / 1000)
    hop_samples = int(args.update_ms * fs / 1000)
    capacity = max(10 * fs, 4 * frame_samples)
    bounds = np.asarray(args.bounds, dtype=float).reshape(3, 2)

    queues = [queue.Queue() for _ in range(3)]
    buffers = [LiveBuffer(channels, capacity) for channels in args.channels]
    streams = []

    for device, channels, q in zip(args.devices, args.channels, queues):
        def callback(indata, frames, time_info, status, q=q):
            if status:
                print(f"Audio status: {status}", flush=True)
            q.put(indata.copy())

        streams.append(
            sd.InputStream(
                device=device,
                channels=channels,
                samplerate=fs,
                dtype="float32",
                blocksize=int(0.02 * fs),
                latency="high",
                callback=callback,
            )
        )

    print("Starting streams. Play equidistant white noise for calibration.")
    for stream in streams:
        stream.start()

    plot = LivePlot(positions, args.channels, bounds) if args.plot else None
    wait = plot.pause if plot else time.sleep

    sync_record = [[] for _ in range(3)]
    deadline = time.monotonic() + args.calibration_seconds
    try:
        while time.monotonic() < deadline:
            drain_queues(queues, buffers, sync_record, args.sync_channels)
            if plot:
                plot.message(
                    f"Calibrating: {deadline - time.monotonic():.0f} s left - "
                    "play loud white noise"
                )
            wait(0.02)
        drain_queues(queues, buffers, sync_record, args.sync_channels)
        if plot:
            plot.message("Computing calibration...")
            plot.pause(0.01)

        # Every stream's sync record starts at its sample 0, the same origin the
        # live buffers use, so the estimated offsets apply to them directly.
        sync_signals = [np.concatenate(blocks) for blocks in sync_record]
        windows = max(6, int(args.calibration_seconds / 2))
        calibrations = [SyncResult(0.0, 0.0, 0.0)]
        for i in (1, 2):
            calibration = estimate_sync(
                sync_signals[0],
                sync_signals[i],
                fs,
                windows=windows,
                window_seconds=0.5,
                max_offset=0.5,
            )
            calibrations.append(calibration)
            print(f"Array {i} calibration: {calibration}", flush=True)
            if calibration.residual_rms > 5e-5:
                print(
                    f"WARNING: array {i} sync is unreliable "
                    f"(residual {calibration.residual_rms * 1e6:.0f} us, want < 20 us). "
                    "Restart and play LOUD white noise during calibration.",
                    flush=True,
                )

        # Leave out channels that carried (almost) no signal during calibration:
        # silent channels only add noise to the SRP-PHAT score.
        levels = np.concatenate(
            [np.sqrt(np.mean(buffer.data**2, axis=0)) for buffer in buffers]
        )
        active = levels > 1e-6
        if not np.all(active):
            dead = [
                f"array {a} ch {c + 1}"
                for a, (lo, n) in enumerate(
                    zip(np.cumsum([0, *args.channels[:-1]]), args.channels)
                )
                for c in range(n)
                if not active[lo + c]
            ]
            print(f"Ignoring {len(dead)} silent channels: {', '.join(dead)}")
        if np.count_nonzero(active) < 4:
            raise RuntimeError("Fewer than 4 working microphones; check the arrays")

        localizer = SRPPhatLocalizer(
            positions[active],
            fs,
            bounds,
            coarse_step=args.coarse_step,
            refinements=args.refinements,
            interpolation=2,
        )

        print("Calibration complete. Localization is running; Ctrl+C to stop.")
        if plot:
            plot.message("Calibration complete - make a sound")
        next_update = 0
        radius = 8

        while True:
            drain_queues(queues, buffers)
            stop = buffers[0].end
            for i in (1, 2):
                scale = 1.0 + calibrations[i].drift
                offset_samples = calibrations[i].offset * fs
                stop = min(
                    stop,
                    int(np.floor((buffers[i].end - radius) * scale - offset_samples)) + 1,
                )
            start = stop - frame_samples

            valid = start >= buffers[0].base and stop >= next_update
            for i in (1, 2):
                scale = 1.0 + calibrations[i].drift
                source_position = (start + calibrations[i].offset * fs) / scale
                valid = valid and source_position >= buffers[i].base + radius - 1

            if valid:
                reference_start = start - buffers[0].base
                reference_stop = stop - buffers[0].base
                synchronized = [buffers[0].data[reference_start:reference_stop]]

                for i in (1, 2):
                    scale = 1.0 + calibrations[i].drift
                    adjusted = SyncResult(
                        calibrations[i].offset - buffers[i].base * scale / fs,
                        calibrations[i].drift,
                        calibrations[i].residual_rms,
                    )
                    synchronized.append(
                        align_array(
                            buffers[i].data,
                            fs,
                            adjusted,
                            frame_samples,
                            start_sample=start,
                        )
                    )

                audio = np.concatenate(synchronized, axis=1)[:, active]
                result = localizer.localize(audio)
                x, y, z = result.position
                level = 20 * np.log10(np.sqrt(np.mean(audio**2)) + 1e-12)
                clear = result.score >= args.min_score
                print(
                    f"Source: x={x:+.3f} m  y={y:+.3f} m  z={z:+.3f} m  "
                    f"score={result.score:.3f}  level={level:5.1f} dBFS"
                    + ("" if clear else "  (no clear source)"),
                    flush=True,
                )
                if plot:
                    plot.update(result.position, result.score, level, clear)
                next_update = stop + hop_samples
            else:
                wait(0.01)

    except KeyboardInterrupt:
        print("\nStopped.")
    finally:
        for stream in streams:
            stream.stop()
            stream.close()


if __name__ == "__main__":
    main()

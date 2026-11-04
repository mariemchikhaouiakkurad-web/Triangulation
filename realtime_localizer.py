"""Real-time three-array SRP-PHAT localization.

Example:
  python realtime_localizer.py --list-devices
  python realtime_localizer.py --devices 1 2 3 --channels 8 16 16 \
      --positions microphone_positions.csv

Keep three_array_localizer.py in the same folder. Play stationary white noise
during calibration. Supply --speaker-position X Y for compensated cross-array
fusion; without it, only board-local pairs are used. Stop with Ctrl+C.
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
    """Unified horizontal source direction and microphone layout."""

    def __init__(self, positions: np.ndarray, channels, bounds: np.ndarray) -> None:
        import matplotlib.pyplot as plt

        self.plt = plt
        plt.ion()
        self.fig, self.top = plt.subplots(figsize=(7, 7))
        names = ["UMA-8", "UMA-16 L", "UMA-16 R"]
        edges = np.cumsum([0, *channels])
        for name, lo, hi in zip(names, edges[:-1], edges[1:]):
            p = positions[lo:hi]
            self.top.scatter(p[:, 0], p[:, 1], s=18, label=name)
        self.top.set_xlabel("x [m] (right)")
        self.top.set_ylabel("y [m] (front)")
        self.top.set_aspect("equal", adjustable="box")
        self.top.grid(alpha=0.3)
        lo = np.minimum(bounds[:2, 0], positions[:, :2].min(axis=0) - 0.1)
        hi = np.maximum(bounds[:2, 1], positions[:, :2].max(axis=0) + 0.1)
        self.top.set_xlim(lo[0], hi[0])
        self.top.set_ylim(lo[1], hi[1])
        self.top.set_title("Unified 2D direction (fixed display radius)")
        self.top.legend(loc="lower right", fontsize=8)
        x_span = np.diff(self.top.get_xlim())[0]
        y_span = np.diff(self.top.get_ylim())[0]
        self.display_radius = 0.35 * min(x_span, y_span)
        angles = np.linspace(0.0, 2.0 * np.pi, 361)
        self.top.plot(
            self.display_radius * np.sin(angles),
            self.display_radius * np.cos(angles),
            "--", color="0.85", lw=1,
        )
        (self.trail,) = self.top.plot([], [], "-", color="0.6", lw=1)
        (self.direction,) = self.top.plot([], [], "-", color="red", lw=2)
        (self.dot,) = self.top.plot([], [], "o", color="red", ms=12)
        self.track: list[np.ndarray] = []
        self.fig.tight_layout()
        self.pause(0.05)

    def update(self, position: np.ndarray | None, score: float, level: float,
               clear: bool) -> None:
        if position is None:
            self.direction.set_data([], [])
            self.dot.set_data([], [])
            self.fig.suptitle(
                f"No accepted direction   score={score:.2f}   "
                f"level={level:.0f} dBFS   no clear sound in selected band"
            )
            return

        x, y = position[:2]
        azimuth = np.degrees(np.arctan2(x, y))
        angle = np.arctan2(x, y)
        display = self.display_radius * np.array([np.sin(angle), np.cos(angle)])
        if clear:
            self.track = (self.track + [display])[-20:]
            points = np.asarray(self.track)
            self.trail.set_data(points[:, 0], points[:, 1])
        self.direction.set_data([0.0, display[0]], [0.0, display[1]])
        self.direction.set_color("red" if clear else "0.6")
        self.dot.set_data([display[0]], [display[1]])
        self.dot.set_color("red" if clear else "0.6")
        status = "live" if clear else "last valid - no clear sound in selected band"
        self.fig.suptitle(
            f"Direction {azimuth:+.2f} deg   score={score:.2f}   "
            f"level={level:.0f} dBFS   {status}"
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


def parse_channels(text: str, available: int) -> list[int]:
    channels: list[int] = []
    for part in text.split(","):
        fields = part.split("-", 1)
        start = int(fields[0])
        stop = int(fields[-1])
        selected = range(start, stop + 1)
        for channel in selected:
            if not 0 <= channel < available:
                raise ValueError(f"Channel {channel} is outside 0-{available - 1}")
            channels.append(channel)
    return sorted(set(channels))


def level_dbfs(
    audio: np.ndarray,
    sample_rate: float,
    frequency_range: tuple[float, float] | None = None,
) -> float:
    if frequency_range is not None:
        low, high = frequency_range
        frequencies = np.fft.rfftfreq(len(audio), d=1.0 / sample_rate)
        band = (frequencies >= low) & (frequencies <= high)
        spectra = np.fft.rfft(audio, axis=0)
        spectra[~band] = 0.0
        audio = np.fft.irfft(spectra, n=len(audio), axis=0)
    rms = np.sqrt(np.mean(audio**2))
    return float(20.0 * np.log10(rms + 1e-12))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list-devices", action="store_true")
    parser.add_argument("--devices", type=int, nargs=3, metavar=("A0", "A1", "A2"))
    parser.add_argument("--channels", type=int, nargs=3, metavar=("N0", "N1", "N2"))
    parser.add_argument("--positions", help="CSV: one x,y,z row per captured channel")
    parser.add_argument("--sync-channels", type=int, nargs=3, default=[0, 0, 0])
    parser.add_argument("--array-channels", nargs=3, default=["0-6", "0-15", "0,14"],
                        metavar=("UMA8", "LEFT", "RIGHT"))
    parser.add_argument("--speaker-position", type=float, nargs=2, metavar=("X", "Y"),
                        help="measured calibration-speaker position; required for cross-array fusion")
    parser.add_argument("--max-sync-error", type=float,
                        help="maximum 95%% synchronization prediction error in seconds (default: one sample)")
    parser.add_argument("--fs", type=int, default=48000)
    parser.add_argument("--calibration-seconds", type=float, default=45.0)
    parser.add_argument("--frame-ms", type=int, default=200)
    parser.add_argument("--update-ms", type=int, default=500)
    parser.add_argument("--coarse-step", type=float, default=0.10)
    parser.add_argument("--refinements", type=int, default=3)
    parser.add_argument("--min-score", type=float, default=0.25,
                        help="estimates below this SRP score are treated as no clear source")
    parser.add_argument("--min-level-db", type=float, default=-50.0,
                        help="minimum level in the selected band for a live result in dBFS (default: -50)")
    parser.add_argument("--frequency-range", type=float, nargs=2,
                        metavar=("LOW_HZ", "HIGH_HZ"),
                        help="use only this frequency band for direction and sound level")
    parser.add_argument("--no-level-gate", action="store_true",
                        help="ignore sound energy; use only the selected band and SRP score")
    parser.add_argument("--plot", action="store_true",
                        help="show the live unified 2D direction and microphone layout")
    parser.add_argument(
        "--bounds", type=float, nargs=6, default=[-1, 1, -1, 1, 0, 0],
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
    if args.frequency_range is not None:
        low, high = args.frequency_range
        if not 0.0 <= low < high <= args.fs / 2:
            parser.error("frequency range must satisfy 0 <= LOW_HZ < HIGH_HZ <= fs/2")

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
    capture_fault = [False] * 3
    streams = []

    for array, (device, channels, q) in enumerate(zip(args.devices, args.channels, queues)):
        def callback(indata, frames, time_info, status, q=q, array=array):
            if status:
                capture_fault[array] = True
                print(f"Audio status on array {array}: {status}", flush=True)
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

    if args.speaker_position is None:
        print("Starting streams. Play white noise for board-local calibration.")
    else:
        print(
            f"Starting streams. Keep the speaker at x={args.speaker_position[0]:.3f}, "
            f"y={args.speaker_position[1]:.3f} m and play white noise."
        )
    for stream in streams:
        stream.start()

    display_edges = np.cumsum([0, *args.channels])
    display_positions = np.vstack([
        positions[lo + np.asarray(parse_channels(selection, count))]
        for lo, selection, count in zip(display_edges[:-1], args.array_channels, args.channels)
    ])
    display_channels = [len(parse_channels(selection, count))
                        for selection, count in zip(args.array_channels, args.channels)]
    plot = LivePlot(display_positions, display_channels, bounds) if args.plot else None
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

        sync_signals = [np.concatenate(blocks) for blocks in sync_record]
        windows = max(6, int(args.calibration_seconds / 2))
        edges = np.cumsum([0, *args.channels])
        sync_indices = edges[:-1] + np.asarray(args.sync_channels)
        sync_positions = positions[sync_indices, :2]
        acoustic_delays = np.zeros(3)
        if args.speaker_position is not None:
            speaker = np.asarray(args.speaker_position)
            distances = np.linalg.norm(sync_positions - speaker, axis=1)
            acoustic_delays = (distances - distances[0]) / 343.0
        calibrations = [SyncResult(0.0, 0.0, 0.0)]
        for i in (1, 2):
            calibration = estimate_sync(
                sync_signals[0],
                sync_signals[i],
                fs,
                windows=windows,
                window_seconds=0.5,
                max_offset=0.5,
                acoustic_delay=acoustic_delays[i],
                max_sync_error=args.max_sync_error,
            )
            calibrations.append(calibration)
            print(f"Array {i} calibration: {calibration}", flush=True)

        use_cross = bool(
            args.speaker_position is not None
            and all(calibration.accepted for calibration in calibrations[1:])
            and not any(capture_fault)
        )
        if args.speaker_position is None:
            print("No measured speaker position: using board-local fusion.")
        elif not use_cross:
            print("Cross-array synchronization rejected: using board-local fusion.")
        else:
            print("Cross-array synchronization accepted.")

        levels = np.concatenate(
            [np.sqrt(np.mean(buffer.data**2, axis=0)) for buffer in buffers]
        )
        selected = np.zeros(total_channels, dtype=bool)
        array_ids = np.zeros(total_channels, dtype=int)
        for array, (lo, hi, count) in enumerate(zip(edges[:-1], edges[1:], args.channels)):
            channels = parse_channels(args.array_channels[array], count)
            selected[lo + np.asarray(channels)] = True
            array_ids[lo:hi] = array
        active = selected & (levels > 1e-6)
        dead = [
            f"array {array} ch {channel + 1}"
            for array, (lo, count) in enumerate(zip(edges[:-1], args.channels))
            for channel in parse_channels(args.array_channels[array], count)
            if not active[lo + channel]
        ]
        if dead:
            print(f"Ignoring silent channels: {', '.join(dead)}")
        if np.count_nonzero(active) < 4:
            raise RuntimeError("Fewer than 4 working microphones; check the arrays")

        def build_localizer(cross_array: bool) -> SRPPhatLocalizer:
            return SRPPhatLocalizer(
                positions[active], fs, bounds,
                coarse_step=args.coarse_step,
                refinements=args.refinements,
                interpolation=2,
                array_ids=array_ids[active],
                use_cross_array=cross_array,
                frequency_range=args.frequency_range,
            )

        localizer = build_localizer(use_cross)

        print("Calibration complete. Localization is running; Ctrl+C to stop.")
        if plot:
            plot.message("Calibration complete - make a sound")
        next_update = 0
        radius = 8
        last_result = None
        last_result_time = 0.0

        while True:
            drain_queues(queues, buffers)
            if any(capture_fault):
                if use_cross:
                    print("Audio discontinuity: switching to board-local fusion.", flush=True)
                    use_cross = False
                    localizer = build_localizer(False)
                capture_fault[:] = [False] * 3
                continue

            synchronized = []
            if use_cross:
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
                    synchronized.append(buffers[0].data[start - buffers[0].base:stop - buffers[0].base])
                    for i in (1, 2):
                        scale = 1.0 + calibrations[i].drift
                        adjusted = SyncResult(
                            calibrations[i].offset - buffers[i].base * scale / fs,
                            calibrations[i].drift,
                            calibrations[i].residual_rms,
                        )
                        synchronized.append(align_array(
                            buffers[i].data, fs, adjusted, frame_samples, start_sample=start
                        ))
            else:
                stop = buffers[0].end
                valid = stop >= next_update
                for buffer in buffers:
                    array_stop = buffer.end
                    array_start = array_stop - frame_samples
                    valid = valid and array_start >= buffer.base
                    if array_start >= buffer.base:
                        synchronized.append(buffer.data[array_start - buffer.base:array_stop - buffer.base])
            if not valid:
                wait(0.01)
                continue

            audio = np.concatenate(synchronized, axis=1)[:, active]
            result = localizer.localize(audio)
            level = level_dbfs(audio, fs, args.frequency_range)
            loud_enough = args.no_level_gate or level >= args.min_level_db
            clear = loud_enough and result.score >= args.min_score
            reason = "no clear sound" if loud_enough else "below level threshold"
            now = time.monotonic()
            if clear:
                last_result = result
                last_result_time = now
                shown = result
                status = "live"
            elif last_result is not None:
                shown = last_result
                status = f"last valid {now - last_result_time:.1f}s ago; {reason}"
            else:
                shown = None
                status = f"no accepted direction; {reason}"
            mode = "cross-array" if use_cross else "board-local"
            direction = f"{shown.azimuth_deg:+7.2f} deg" if shown is not None else "--"
            print(
                f"Direction: {direction}  "
                f"score={result.score:.3f} level={level:5.1f} dBFS  "
                f"{mode}; {status}",
                flush=True,
            )
            if plot:
                plot.update(
                    shown.position if shown is not None else None,
                    result.score, level, clear,
                )
            next_update = stop + hop_samples

    except KeyboardInterrupt:
        print("\nStopped.")
    finally:
        for stream in streams:
            stream.stop()
            stream.close()


if __name__ == "__main__":
    main()

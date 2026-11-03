"""Synchronize three microphone arrays, locate a source, and beamform it.

Algorithms
----------
1. Windowed GCC-PHAT + linear regression: capture offset/drift estimation.
2. Fractional-delay resampling: alignment to the reference array's timeline.
3. Coarse-to-fine SRP-PHAT: robust 2-D/3-D source localization.
4. Frequency-domain fractional-delay Delay-and-Sum: source enhancement.

Important
---------
* Best synchronization: all arrays share one hardware sample clock. If they are
  independent USB devices, use a common electrical calibration pulse/loopback.
  Acoustic calibration also contains propagation delays and can bias the result.
* Microphone coordinates must be in one right-handed coordinate system, in
  metres, and in exactly the same order as the concatenated audio channels.
* A planar layout cannot distinguish mirrored sources above/below its plane.
  Restrict search_bounds to the physically valid half-space.

The routines require NumPy and SciPy. Audio arrays use shape (samples, channels).
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Sequence

import numpy as np
from numpy.typing import NDArray
from scipy import fft

FloatArray = NDArray[np.floating]


@dataclass
class SyncResult:
    """Timing of one array relative to the reference array.

    offset is in seconds; positive means the array initially lags the reference.
    drift is the fractional sample-time scale error (1e-6 equals 1 ppm).
    residual_rms measures calibration quality and should be near zero for a
    valid electrical/common-clock calibration.
    """

    offset: float
    drift: float
    residual_rms: float


@dataclass
class LocalizationResult:
    position: NDArray[np.float64]
    score: float


def gcc_phat(
    signal: NDArray,
    reference: NDArray,
    sample_rate: float,
    max_delay: float,
) -> tuple[float, float]:
    """Return (delay_seconds, peak) using interpolated GCC-PHAT.

    A positive delay means that signal lags reference.
    """
    signal = np.asarray(signal, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    if signal.ndim != 1 or reference.ndim != 1 or len(signal) != len(reference):
        raise ValueError("GCC-PHAT inputs must be equally sized 1-D signals")

    n = len(signal)
    nfft = fft.next_fast_len(2 * n)
    window = np.hanning(n)
    spectrum_s = fft.rfft(signal * window, nfft)
    spectrum_r = fft.rfft(reference * window, nfft)
    cross = spectrum_s * np.conj(spectrum_r)
    cross /= np.maximum(np.abs(cross), 1e-15)

    correlation = fft.irfft(cross, nfft)
    max_shift = min(int(np.ceil(max_delay * sample_rate)), n // 2 - 1)
    shifts = np.arange(-max_shift, max_shift + 1)
    values = correlation[shifts % nfft]
    peak_index = int(np.argmax(values))

    # Three-point parabolic interpolation gives a fractional-sample estimate.
    fractional = 0.0
    if 0 < peak_index < len(values) - 1:
        left, centre, right = values[peak_index - 1 : peak_index + 2]
        denominator = left - 2.0 * centre + right
        if abs(denominator) > 1e-20:
            fractional = 0.5 * (left - right) / denominator
            fractional = float(np.clip(fractional, -1.0, 1.0))

    delay = (shifts[peak_index] + fractional) / sample_rate
    return float(delay), float(values[peak_index])


def estimate_sync(
    reference: NDArray,
    other: NDArray,
    sample_rate: float,
    *,
    windows: int = 8,
    window_seconds: float = 0.5,
    max_offset: float = 0.25,
) -> SyncResult:
    """Estimate constant offset and linear clock drift from a calibration recording.

    reference and other are single channels containing the same electrical pulse
    or common calibration stimulus. At least two separated windows are required;
    more windows make the drift estimate more reliable.
    """
    reference = np.asarray(reference, dtype=np.float64).squeeze()
    other = np.asarray(other, dtype=np.float64).squeeze()
    available = min(len(reference), len(other))
    window_length = int(window_seconds * sample_rate)
    if available < 2 * window_length or windows < 2:
        raise ValueError("Recording is too short for offset/drift estimation")

    starts = np.linspace(0, available - window_length, windows, dtype=int)
    times: list[float] = []
    delays: list[float] = []
    for start in starts:
        stop = start + window_length
        delay, _ = gcc_phat(
            other[start:stop],
            reference[start:stop],
            sample_rate,
            max_offset,
        )
        times.append((start + 0.5 * window_length) / sample_rate)
        delays.append(delay)

    slope, offset = np.polyfit(times, delays, 1)
    fitted = offset + slope * np.asarray(times)
    residual = float(np.sqrt(np.mean((np.asarray(delays) - fitted) ** 2)))

    # If measured delay(t) = offset - drift*t, the array's sample-time scale
    # is (1 + drift). Resampling below uses that same convention.
    return SyncResult(float(offset), float(-slope), residual)


def align_array(
    audio: NDArray,
    sample_rate: float,
    calibration: SyncResult,
    output_samples: int,
    *,
    start_sample: int = 0,
) -> FloatArray:
    """Resample an array onto the reference timeline.

    Chunked Lanczos-8 interpolation provides band-limited fractional delays and
    supports the slowly changing sample position caused by clock drift.
    """
    audio = np.asarray(audio, dtype=np.float64)
    if audio.ndim == 1:
        audio = audio[:, None]
    if audio.ndim != 2:
        raise ValueError("audio must have shape (samples, channels)")

    output_time = np.arange(
        start_sample, start_sample + output_samples, dtype=np.float64
    )
    source_positions = (
        output_time + calibration.offset * sample_rate
    ) / (1.0 + calibration.drift)

    kernel_radius = 8
    if (
        source_positions[0] < kernel_radius - 1
        or source_positions[-1] > len(audio) - kernel_radius
    ):
        raise ValueError(
            "Calibration requires samples outside the recording; trim the "
            "reference output or capture additional guard samples"
        )

    aligned = np.empty((output_samples, audio.shape[1]), dtype=np.float64)
    taps = np.arange(-kernel_radius + 1, kernel_radius + 1)
    chunk_size = 65_536
    for start in range(0, output_samples, chunk_size):
        stop = min(start + chunk_size, output_samples)
        positions = source_positions[start:stop]
        indices = np.floor(positions).astype(np.int64)[:, None] + taps[None, :]
        distances = positions[:, None] - indices

        # Lanczos interpolation: a sinc delayed by a fractional sample, windowed
        # by a second sinc. Normalization preserves unity gain at DC.
        weights = np.sinc(distances) * np.sinc(distances / kernel_radius)
        weights /= np.sum(weights, axis=1, keepdims=True)
        for channel in range(audio.shape[1]):
            aligned[start:stop, channel] = np.sum(
                audio[indices, channel] * weights, axis=1
            )

    return aligned


def synchronize_arrays(
    arrays: Sequence[NDArray],
    sample_rate: float,
    sync_channels: Sequence[int],
    *,
    windows: int = 8,
    window_seconds: float = 0.5,
    max_offset: float = 0.25,
) -> tuple[FloatArray, list[SyncResult]]:
    """Synchronize arrays to arrays[0], returning one concatenated audio matrix.

    Use recordings of a common electrical calibration event. Each entry of
    sync_channels selects the channel carrying that event in its array.
    """
    if len(arrays) != 3 or len(sync_channels) != 3:
        raise ValueError("This configuration expects exactly three arrays")

    reference = np.asarray(arrays[0], dtype=np.float64)
    other_arrays = [np.asarray(array, dtype=np.float64) for array in arrays[1:]]
    results = [SyncResult(0.0, 0.0, 0.0)]
    reference_sync = reference[:, sync_channels[0]]

    for array, channel in zip(other_arrays, sync_channels[1:]):
        results.append(
            estimate_sync(
                reference_sync,
                array[:, channel],
                sample_rate,
                windows=windows,
                window_seconds=window_seconds,
                max_offset=max_offset,
            )
        )

    # Find the interval for which every resampler has all required input taps.
    # This automatically removes invalid edges caused by offset and drift.
    radius = 8
    start = 0
    stop = len(reference)
    for array, calibration in zip(other_arrays, results[1:]):
        scale = 1.0 + calibration.drift
        offset_samples = calibration.offset * sample_rate
        start = max(start, int(np.ceil((radius - 1) * scale - offset_samples)))
        stop = min(
            stop,
            int(np.floor((len(array) - radius) * scale - offset_samples)) + 1,
        )
    if stop - start < 1024:
        raise ValueError("Less than 1024 common valid samples remain after alignment")

    synchronized = [reference[start:stop]]
    for array, calibration in zip(other_arrays, results[1:]):
        synchronized.append(
            align_array(
                array,
                sample_rate,
                calibration,
                stop - start,
                start_sample=start,
            )
        )

    return np.concatenate(synchronized, axis=1), results


class SRPPhatLocalizer:
    """Coarse-to-fine Steered-Response-Power PHAT source localizer."""

    def __init__(
        self,
        microphone_positions: NDArray,
        sample_rate: float,
        search_bounds: NDArray,
        *,
        speed_of_sound: float = 343.0,
        coarse_step: float = 0.05,
        refinements: int = 3,
        interpolation: int = 4,
    ) -> None:
        self.positions = np.asarray(microphone_positions, dtype=np.float64)
        self.sample_rate = float(sample_rate)
        self.bounds = np.asarray(search_bounds, dtype=np.float64)
        self.c = float(speed_of_sound)
        self.coarse_step = float(coarse_step)
        self.refinements = int(refinements)
        self.interpolation = int(interpolation)

        if self.positions.ndim != 2 or self.positions.shape[1] != 3:
            raise ValueError("microphone_positions must have shape (microphones, 3)")
        if self.bounds.shape != (3, 2) or np.any(self.bounds[:, 0] > self.bounds[:, 1]):
            raise ValueError("search_bounds must be [[xmin,xmax],[ymin,ymax],[zmin,zmax]]")

        self.pairs = np.asarray(list(combinations(range(len(self.positions)), 2)))
        diameters = np.linalg.norm(
            self.positions[:, None, :] - self.positions[None, :, :], axis=-1
        )
        self.max_tdoa = float(np.max(diameters) / self.c + 2.0 / self.sample_rate)
        self._delay_grid: NDArray | None = None
        self._pair_responses: NDArray | None = None

    def _prepare_responses(self, audio: NDArray) -> None:
        """Precompute interpolated pairwise GCC-PHAT functions for one frame."""
        audio = np.asarray(audio, dtype=np.float64)
        if audio.ndim != 2 or audio.shape[1] != len(self.positions):
            raise ValueError("Audio channels do not match microphone coordinates")

        audio = audio - np.mean(audio, axis=0, keepdims=True)
        audio *= np.hanning(len(audio))[:, None]
        nfft = fft.next_fast_len(2 * len(audio))
        spectra = fft.rfft(audio, nfft, axis=0)

        upsampled_size = nfft * self.interpolation
        max_shift = int(np.ceil(self.max_tdoa * self.sample_rate * self.interpolation))
        shifts = np.arange(-max_shift, max_shift + 1)
        self._delay_grid = shifts / (self.sample_rate * self.interpolation)
        self._pair_responses = np.empty(
            (len(self.pairs), len(shifts)), dtype=np.float32
        )

        for pair_number, (first, second) in enumerate(self.pairs):
            cross = spectra[:, first] * np.conj(spectra[:, second])
            cross /= np.maximum(np.abs(cross), 1e-15)
            correlation = fft.irfft(cross, upsampled_size)
            response = correlation[shifts % upsampled_size]
            response /= max(float(np.max(np.abs(response))), 1e-15)
            self._pair_responses[pair_number] = response

    @staticmethod
    def _grid(bounds: NDArray, step: float) -> NDArray:
        axes = []
        for lower, upper in bounds:
            if np.isclose(lower, upper):
                axes.append(np.asarray([lower]))
            else:
                axes.append(np.arange(lower, upper + 0.5 * step, step))
        mesh = np.meshgrid(*axes, indexing="ij")
        return np.column_stack([axis.ravel() for axis in mesh])

    def _score(self, candidates: NDArray) -> NDArray:
        assert self._delay_grid is not None and self._pair_responses is not None
        distances = np.linalg.norm(
            candidates[:, None, :] - self.positions[None, :, :], axis=-1
        )
        predicted_tdoa = (
            distances[:, self.pairs[:, 0]] - distances[:, self.pairs[:, 1]]
        ) / self.c

        scores = np.zeros(len(candidates), dtype=np.float64)
        for pair_number in range(len(self.pairs)):
            scores += np.interp(
                predicted_tdoa[:, pair_number],
                self._delay_grid,
                self._pair_responses[pair_number],
                left=0.0,
                right=0.0,
            )
        return scores / len(self.pairs)

    def localize(self, audio_frame: NDArray) -> LocalizationResult:
        """Locate the dominant broadband source in a short synchronized frame.

        A frame of roughly 0.1-0.3 seconds is a useful starting point. Track
        speech by applying this method to overlapping frames and optionally
        restricting subsequent searches around the previous valid estimate.
        """
        self._prepare_responses(audio_frame)
        bounds = self.bounds.copy()
        step = self.coarse_step
        best_position = np.mean(bounds, axis=1)
        best_score = -np.inf

        for level in range(self.refinements + 1):
            candidates = self._grid(bounds, step)
            if len(candidates) > 500_000:
                raise ValueError(
                    f"Search grid has {len(candidates):,} points; increase "
                    "coarse_step or reduce search_bounds"
                )
            scores = self._score(candidates)
            best = int(np.argmax(scores))
            best_position = candidates[best]
            best_score = float(scores[best])

            if level < self.refinements:
                radius = 1.5 * step
                bounds = np.column_stack(
                    (
                        np.maximum(self.bounds[:, 0], best_position - radius),
                        np.minimum(self.bounds[:, 1], best_position + radius),
                    )
                )
                step /= 4.0

        return LocalizationResult(np.asarray(best_position), best_score)


def delay_and_sum(
    audio: NDArray,
    microphone_positions: NDArray,
    source_position: NDArray,
    sample_rate: float,
    *,
    speed_of_sound: float = 343.0,
) -> FloatArray:
    """Delay-and-Sum beamform synchronized audio toward source_position."""
    audio = np.asarray(audio, dtype=np.float64)
    positions = np.asarray(microphone_positions, dtype=np.float64)
    source_position = np.asarray(source_position, dtype=np.float64)
    if audio.shape[1] != len(positions):
        raise ValueError("Audio channels do not match microphone coordinates")

    propagation_delays = np.linalg.norm(positions - source_position, axis=1) / speed_of_sound
    # Remove the common delay; only relative delays matter for beamforming.
    propagation_delays -= np.min(propagation_delays)

    nfft = fft.next_fast_len(len(audio))
    frequencies = fft.rfftfreq(nfft, d=1.0 / sample_rate)
    spectra = fft.rfft(audio, nfft, axis=0)
    steering = np.exp(2j * np.pi * frequencies[:, None] * propagation_delays[None, :])
    aligned = spectra * steering
    output = fft.irfft(np.mean(aligned, axis=1), nfft)[: len(audio)]
    return np.asarray(output, dtype=np.float64)


# Minimal integration example (replace the placeholders with measured data):
#
# fs = 48_000
# raw_arrays = [array0, array1, array2]       # each: samples x channels
# microphones = np.array([...])               # total channels x 3, metres
# bounds = np.array([[-1, 1], [-1, 1], [0.05, 2.0]])
#
# audio, sync = synchronize_arrays(raw_arrays, fs, sync_channels=[0, 16, 32])
# localizer = SRPPhatLocalizer(microphones, fs, bounds, coarse_step=0.05)
# estimate = localizer.localize(audio[: int(0.20 * fs)])
# print("Estimated source [m]:", estimate.position, "SRP score:", estimate.score)
# enhanced = delay_and_sum(audio, microphones, estimate.position, fs)

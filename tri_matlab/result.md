# Frequency mode (`freq_live`) – analysis of session 20261102_163609

## 1. The drone sound is mostly not detected, so no direction is computed
- Only **50 of 1476 updates (3 %)** detected the harmonic comb. A direction is only computed
  after a detection; the spectrum always updates, which is why the frequency is seen moving
  but no direction appears.
- **Cause:** the test that rejects false alarms (room / phone-speaker resonances) requires every
  harmonic to keep a **steady frequency over the 1 s window** (phase coherence ≥ 0.85). A drone
  sound that **changes frequency** fails it:
  - t = 120 s: several peaks at coherence 0.33–0.65
  - t = 200 s: the strongest tone (3428 Hz, 67 dB above noise) only 0.34
- The test protects well against false alarms but is too strict for a drone whose speed varies.

## 2. When it detects, the fundamental is wrong (bug)
- At 75–89 s there is a clean comb **1137, 2273, 3410, (5684) Hz** → **f0 = 1137 Hz**,
  narrow and steady.
- The detector reported **f0 = 103.2 Hz**: 1137 Hz is the 11th harmonic of 103.2 Hz, so that
  candidate explains the same peaks equally well. The octave check only handles 2× and 3×.
- With the wrong f0, low-frequency bins (< 400 Hz) are used. There the wavelength (0.9–3.4 m) is
  much larger than the boards (13–18 cm), hence the very wide directions (**±42° to ±182°**).

## 3. To check: does the 1137 Hz family really come from the phone?
- At 75–89 s the overall level was very low (−64 to −67 dBFS) but this tone was 47–49 dB above the
  noise and perfectly steady. It may be the phone at low volume, or **electronic whine** (USB
  board, power supply, laptop), which has no meaningful direction.
- **Test:** record ~20 s with the phone silent. If 1137 Hz is still there, it is electrical.

## Calibration in the frequency mode
- The calibration was not removed: `freq_live` **reuses the calibration saved by the TDOA mode**
  (`tdoa_calibration.mat`, board orientations of 2026-11-02 15:40). Orientations only change when a
  board is moved.
- **If a board was moved since 15:40, the saved orientations are wrong and every direction is
  wrong.** Recalibrate with `cfg.mode = 'tdoa_live'`, `cfg.tdoa.reuse_cal = false`.

## Proposed fixes (not applied yet)
1. **Correct fundamental:** prefer the highest f0 that explains the peaks (its first harmonics must
   be present), so 1137 Hz is not reported as 103 Hz.
2. **Accept a drone whose frequency varies:** measure steadiness over a shorter time and let each
   harmonic drift in proportion to its order; resonances still fail.
3. **Optional calibration at the start of the frequency mode** (`cfg.freq.calibrate = true`).
4. **Ignore harmonics below ~800 Hz for the direction** (boards too small to resolve direction
   there); keep them for detection only.

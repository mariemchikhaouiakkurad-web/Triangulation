#!/usr/bin/env python3
"""
uma16_tap_map.py - guided tap test: which physical mic drives which channel (cut UMA-16 v2).

Prompts you to tap each remaining mic in turn (3 s each, after 3 s of silence), then reports:
  - per mic: the channel(s) that react and by how many dB (real mic: > +15 dB)
  - channel groups carrying the SAME signal (copies / shorted / floating lines)

Usage: python3 uma16_tap_map.py [--device hw:UMA16v2,0] [--file existing.wav]
"""
import argparse
import os
import subprocess
import sys
import time
import wave

import numpy as np

MICS = ['MIC1', 'MIC3', 'MIC5', 'MIC7', 'MIC10', 'MIC12', 'MIC14', 'MIC16']
NCH = 16
SLOT = 3.0          # seconds per step
FS = 48000


def record(device, path, seconds):
    p = subprocess.Popen(['arecord', '-q', '-D', device, '-c', str(NCH), '-f', 'S32_LE', '-r', str(FS),
                          '-d', str(int(seconds) + 1), path])
    t0 = time.time()
    steps = ['SILENCE (do not touch)'] + [f'TAP {m} repeatedly' for m in MICS]
    for i, s in enumerate(steps):
        while time.time() - t0 < i * SLOT:
            time.sleep(0.01)
        print(f'[{i * SLOT:4.0f}s] >>> {s}', flush=True)
    p.wait()
    if p.returncode != 0:
        sys.exit('arecord failed (device busy? close ODAS / other recorders)')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--device', default='hw:UMA16v2,0')
    ap.add_argument('--file', default=None, help='analyse an existing recording made with this script')
    ap.add_argument('--mics', nargs='+', default=None, help='mic labels to tap, in order (e.g. MIC8 MIC9 MIC2 MIC15)')
    ap.add_argument('--nch', type=int, default=16, help='channels of the device (UMA-8 RAW: 8)')
    a = ap.parse_args()
    global MICS, NCH
    if a.mics:
        MICS = a.mics
    NCH = a.nch
    path = a.file or os.path.expanduser('~/uma16_tapmap.wav')
    n_steps = 1 + len(MICS)
    if not a.file:
        print('Get ready: first 3 s silence, then tap each mic when its name appears.')
        time.sleep(2)
        record(a.device, path, n_steps * SLOT)

    w = wave.open(path)
    x = np.frombuffer(w.readframes(w.getnframes()), np.int32).reshape(-1, NCH).astype(float)
    x -= np.median(x, axis=0)
    seg = int(SLOT * FS)
    guard = int(0.3 * FS)                    # ignore 0.3 s at each step boundary (reaction time)

    def level(i):
        s = x[i * seg + guard:(i + 1) * seg - guard]
        return 20 * np.log10(np.sqrt((s ** 2).mean(0)) + 1)
    base = level(0)

    print('\n=== which channel reacts to each tap (dB above silence) ===')
    print('mic     best ch   +dB    other channels > +10 dB')
    for i, m in enumerate(MICS, start=1):
        d = level(i) - base
        best = int(np.argmax(d))
        others = [f'ch{c}({d[c]:+.0f})' for c in np.argsort(-d) if c != best and d[c] > 10][:6]
        verdict = 'OK' if d[best] > 15 else 'NO RESPONSE'
        print(f'{m:6s}  ch{best:<3d}  {d[best]:+5.1f}  {verdict:12s} MATLAB channel {best + 1}   {" ".join(others)}')

    print('\n=== channels carrying the same signal (|corr| > 0.999) ===')
    with np.errstate(invalid='ignore', divide='ignore'):
        c = np.nan_to_num(np.corrcoef(x.T))
    seen = set()
    for i in range(NCH):
        if i in seen:
            continue
        grp = [j for j in range(NCH) if abs(c[i, j]) > 0.999]
        if len(grp) > 1:
            seen.update(grp)
            print('  same signal: ' + ' '.join(f'ch{j}{"(-)" if c[i, j] < 0 else ""}' for j in grp))
    print('\nReading:\n'
          '  pair (2k, 2k+1) identical  = normal here: partner mic was cut off, slot holds a copy\n'
          '  >2 channels identical      = lines shorted together or floating (no mic driving them)\n'
          '  mic with NO RESPONSE       = no clock / no 3.3 V / broken data trace to that mic')


if __name__ == '__main__':
    main()

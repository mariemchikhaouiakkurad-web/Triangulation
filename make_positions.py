"""Write microphone_positions.csv for one UMA-8 in front and two UMA-16.

Coordinate system (metres): origin = UMA-8 centre, +x to the right,
+y forward (towards where the sound is), +z up.

        y (front, towards the sound)
        ^
        |      [UMA-8]  (0, 0)
        |
  [UMA-16 L]          [UMA-16 R]     (-spread, -back)  (+spread, -back)
        +----------------------> x

Example:
  python make_positions.py --spread 0.6 --back 0.3 --facing up

--facing up    arrays lie flat on a table (mics in the x-y plane)
--facing front arrays stand upright facing +y (mics in the x-z plane)

Row order matches --devices: UMA-8, left UMA-16, right UMA-16.

VERIFY the per-array templates below against the miniDSP documentation for
your firmware: the radius, which channel is the centre mic, what the 8th UMA-8
channel carries, and the UMA-16 channel-to-mic order (assumed row-major here).
A wrong channel order gives confident but wrong positions. Test by tapping
single mics and checking which channel responds.
"""

import argparse

import numpy as np

# UMA-8: channel 0 = centre, channels 1-6 = ring, channel 7 placed at centre
# (assumption: check what your firmware outputs on the 8th channel).
ring = np.arange(6) * 2 * np.pi / 6
UMA8 = np.vstack([
    [0, 0, 0],
    np.c_[0.043 * np.cos(ring), 0.043 * np.sin(ring), np.zeros(6)],
    [0, 0, 0],
])

# UMA-16: 4x4 grid, 42 mm pitch, assumed row-major channel order.
g = (np.arange(4) - 1.5) * 0.042
UMA16 = np.array([[x, y, 0.0] for y in g for x in g])


def place(template, center, facing):
    if facing == "front":
        # Upright board: the template's y axis becomes height.
        template = template[:, [0, 2, 1]]
    return template + center


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--spread", type=float, default=0.6,
                        help="left/right distance of each UMA-16 from the UMA-8 (m)")
    parser.add_argument("--back", type=float, default=0.3,
                        help="how far behind the UMA-8 the UMA-16s are (m)")
    parser.add_argument("--height16", type=float, default=0.0,
                        help="height of the UMA-16s relative to the UMA-8 (m)")
    parser.add_argument("--facing", choices=("up", "front"), default="up")
    parser.add_argument("--out", default="microphone_positions.csv")
    args = parser.parse_args()

    centers = [
        [0.0, 0.0, 0.0],                                # UMA-8, front
        [-args.spread, -args.back, args.height16],      # UMA-16, left
        [args.spread, -args.back, args.height16],       # UMA-16, right
    ]
    positions = np.vstack([
        place(t, c, args.facing) for t, c in zip((UMA8, UMA16, UMA16), centers)
    ])
    np.savetxt(args.out, positions, delimiter=",", fmt="%.4f",
               header="x,y,z", comments="")
    print(f"Wrote {len(positions)} rows to {args.out} "
          f"(UMA-8 front, UMA-16 left/right at x=±{args.spread}, y=-{args.back}, "
          f"facing {args.facing})")


if __name__ == "__main__":
    main()

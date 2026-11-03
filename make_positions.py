"""Write microphone_positions.csv for one UMA-8 in front and two UMA-16.

Coordinate system (metres): origin = UMA-8 centre, +x to the right,
+y forward (towards where the sound is), +z up.

        y (front, towards the sound)
        ^
        |      [UMA-8]  (0, 0)
        |
  [UMA-16 L]          [UMA-16 R]     (-spread, -back)  (+spread, -back)
        +----------------------> x

The default photographed layout uses the centers and approximate board angles
documented in tri_matlab/tri_config.m. Use --layout custom for another layout.
Vertical coordinates are removed by default for direction-only operation.

Example:
  python make_positions.py --layout photographed --out microphone_positions.csv

Row order matches --devices: UMA-8, left UMA-16, right UMA-16. The microphone
templates match tri_matlab/uma8.cfg and uma16.cfg. The real-time default uses
UMA-8 channels 0-6, all left channels, and right channels 0 and 14 (MIC1/MIC16).
Verify the channel mapping and board angles on the physical hardware.
"""

import argparse

import numpy as np

UMA8 = np.array([
    [0.0, 0.0, 0.0],
    [-0.02819, 0.01026, 0.0],
    [-0.02298, -0.01928, 0.0],
    [0.00521, -0.02954, 0.0],
    [0.02819, -0.01026, 0.0],
    [0.02298, 0.01928, 0.0],
    [-0.00521, 0.02954, 0.0],
    [0.0, 0.0, 0.0],
])

UMA16 = np.array([
    [-0.021, -0.063, 0.0], [-0.063, -0.063, 0.0],
    [-0.021, -0.021, 0.0], [-0.063, -0.021, 0.0],
    [-0.021, 0.021, 0.0], [-0.063, 0.021, 0.0],
    [-0.021, 0.063, 0.0], [-0.063, 0.063, 0.0],
    [0.063, 0.063, 0.0], [0.021, 0.063, 0.0],
    [0.063, 0.021, 0.0], [0.021, 0.021, 0.0],
    [0.063, -0.021, 0.0], [0.021, -0.021, 0.0],
    [0.063, -0.063, 0.0], [0.021, -0.063, 0.0],
])


def place(template, center, facing, yaw):
    if facing == "front":
        template = template[:, [0, 2, 1]]
    angle = np.deg2rad(yaw)
    rotation = np.array([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])
    return template @ rotation.T + center


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--layout", choices=("photographed", "custom"), default="photographed")
    parser.add_argument("--spread", type=float, default=0.6,
                        help="custom left/right UMA-16 distance from UMA-8 (m)")
    parser.add_argument("--back", type=float, default=0.3,
                        help="custom UMA-16 distance behind UMA-8 (m)")
    parser.add_argument("--height16", type=float, default=0.0,
                        help="custom UMA-16 height relative to UMA-8 (m)")
    parser.add_argument("--facing", choices=("up", "front"), default="front")
    parser.add_argument("--yaw", type=float, nargs=3, default=None,
                        help="custom UMA-8/left/right rotation in the horizontal plane (degrees)")
    parser.add_argument("--keep-height", action="store_true",
                        help="retain vertical coordinates instead of producing a 2D layout")
    parser.add_argument("--out", default="microphone_positions.csv")
    args = parser.parse_args()

    if args.layout == "photographed":
        centers = [[0.0, 0.0, 0.0], [-0.073, -0.239, 0.0], [0.076, -0.249, 0.0]]
        yaw = args.yaw if args.yaw is not None else [0.0, 73.0, -73.0]
    else:
        centers = [
            [0.0, 0.0, 0.0],
            [-args.spread, -args.back, args.height16],
            [args.spread, -args.back, args.height16],
        ]
        yaw = args.yaw if args.yaw is not None else [0.0, 0.0, 0.0]
    positions = np.vstack([
        place(t, c, args.facing, angle)
        for t, c, angle in zip((UMA8, UMA16, UMA16), centers, yaw)
    ])
    if not args.keep_height:
        positions[:, 2] = 0.0
    np.savetxt(args.out, positions, delimiter=",", fmt="%.5f",
               header="x,y,z", comments="")
    print(f"Wrote {len(positions)} rows to {args.out} "
          f"({args.layout} layout, yaw={yaw}, facing={args.facing}, 2D={not args.keep_height})")


if __name__ == "__main__":
    main()

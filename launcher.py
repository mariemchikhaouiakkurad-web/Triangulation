#!/usr/bin/env python3
"""An easy-to-use launcher for the Triangulation project."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
DEFAULT_PROJECT = HERE if (HERE / "realtime_localizer.py").is_file() else HERE.parent / "Triangulation"
DEFAULT_POSITIONS = HERE / "microphone_positions.csv"


def project_from(args: argparse.Namespace) -> Path:
    project = args.project.expanduser().resolve()
    required = (
        "realtime_localizer.py",
        "three_array_localizer.py",
        "make_positions.py",
    )
    missing = [name for name in required if not (project / name).is_file()]
    if missing:
        raise ValueError(f"{project} is not a Triangulation project; missing: {', '.join(missing)}")
    if getattr(args, "positions", None) is None and DEFAULT_POSITIONS.is_file():
        args.positions = DEFAULT_POSITIONS
    return project


def print_command(command: list[str]) -> None:
    print("\nCommand:", shlex.join(command))


def execute(command: list[str], project: Path, dry_run: bool = False) -> int:
    print_command(command)
    if dry_run:
        print("Dry run only; nothing was started.")
        return 0
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    completed = subprocess.run(command, cwd=project, env=env)
    return completed.returncode


def check_command(args: argparse.Namespace) -> int:
    project = project_from(args)
    positions = args.positions.expanduser().resolve() if args.positions else project / "microphone_positions.csv"
    ok = True
    print(f"Project: {project}")
    print("Localization algorithms are invoked without source edits.")
    for package in ("numpy", "scipy", "sounddevice", "matplotlib"):
        installed = importlib.util.find_spec(package) is not None
        print(f"  {package:12} {'OK' if installed else 'MISSING'}")
        ok = ok and installed
    if not positions.is_file():
        print(f"Microphone positions: MISSING ({positions})")
        return 1
    with positions.open(newline="") as handle:
        rows = []
        for row in csv.reader(handle):
            if not row:
                continue
            try:
                values = [float(value) for value in row]
            except ValueError:
                continue
            rows.append(values)
    if len(rows) != 40 or any(len(row) != 3 for row in rows):
        print(f"Microphone positions: INVALID (expected 40 rows of x,y,z; found {len(rows)})")
        ok = False
    else:
        print(f"Microphone positions: OK ({positions})")
        print("  40 captured channels; real-time selection is 7 + 16 + 2 microphones.")
    print("\nSetup is ready." if ok else "\nFix the missing items above before running localization.")
    return 0 if ok else 1


def devices_command(args: argparse.Namespace) -> int:
    project = project_from(args)
    return execute(
        [sys.executable, str(project / "realtime_localizer.py"), "--list-devices"],
        project,
        args.dry_run,
    )


def geometry_command(args: argparse.Namespace) -> int:
    project = project_from(args)
    output = args.out.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        str(project / "make_positions.py"),
        "--layout", args.layout,
        "--spread", str(args.spread),
        "--back", str(args.back),
        "--height16", str(args.height16),
        "--facing", args.facing,
        "--out", str(output),
    ]
    result = execute(command, project, args.dry_run)
    if result == 0 and not args.dry_run:
        print(f"\nUse this layout when starting localization:\n  --positions {output}")
    return result


def prompt(text: str, default: str | None = None) -> str:
    suffix = f" [{default}]" if default is not None else ""
    value = input(f"{text}{suffix}: ").strip()
    return value or (default or "")


def prompt_devices() -> list[int]:
    while True:
        value = prompt("Device IDs for UMA-8, UMA-16 left, UMA-16 right (example: 1 2 3)")
        try:
            devices = [int(part) for part in value.split()]
            if len(devices) != 3:
                raise ValueError
            return devices
        except ValueError:
            print("Enter exactly three integer device IDs. Use 'List audio devices' first if needed.")


def run_command(args: argparse.Namespace) -> int:
    project = project_from(args)
    interactive = sys.stdin.isatty()
    if args.devices is None:
        if not interactive:
            raise ValueError("run requires --devices A0 A1 A2 when input is not interactive")
        args.devices = prompt_devices()
    if args.speaker_position is None and interactive:
        value = prompt("Speaker x y during calibration in metres", "0 1").split()
        if len(value) != 2:
            raise ValueError("speaker position requires x and y")
        args.speaker_position = [float(value[0]), float(value[1])]
    positions = args.positions.expanduser().resolve() if args.positions else project / "microphone_positions.csv"
    if not positions.is_file():
        raise ValueError(f"microphone positions not found: {positions}")
    bounds = list(args.bounds)
    bounds[4] = 0.0
    bounds[5] = 0.0
    command = [
        sys.executable,
        str(project / "realtime_localizer.py"),
        "--devices", *[str(value) for value in args.devices],
        "--channels", *[str(value) for value in args.channels],
        "--positions", str(positions),
        "--sync-channels", *[str(value) for value in args.sync_channels],
        "--array-channels", "0-6", "0-15", "0,14",
        "--fs", str(args.fs),
        "--calibration-seconds", str(args.calibration_seconds),
        "--frame-ms", str(args.frame_ms),
        "--update-ms", str(args.update_ms),
        "--coarse-step", str(args.coarse_step),
        "--refinements", str(args.refinements),
        "--min-score", str(args.min_score),
        "--bounds", *[str(value) for value in bounds],
    ]
    if args.speaker_position is not None:
        command.extend(["--speaker-position", *[str(value) for value in args.speaker_position]])
    if args.max_sync_error is not None:
        command.extend(["--max-sync-error", str(args.max_sync_error)])
    if args.plot:
        command.append("--plot")
    print("\nMode: unified 2D direction (front=0 deg, right=+90 deg, left=-90 deg)")
    print(f"Positions: {positions}")
    if args.speaker_position is not None:
        print(f"Calibration speaker: {args.speaker_position} m")
    print(f"Calibration: play loud white noise for {args.calibration_seconds:g} seconds.")
    print("Stop localization with Ctrl+C.")
    return execute(command, project, args.dry_run)


def add_dry_run(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--dry-run", action="store_true", help="print the command without running it")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog="Run without a command to open the interactive menu.",
    )
    parser.add_argument("--project", type=Path, default=DEFAULT_PROJECT, help=argparse.SUPPRESS)
    subparsers = parser.add_subparsers(dest="command")

    check = subparsers.add_parser("check", help="check Python packages and microphone coordinates")
    check.add_argument("--positions", type=Path)
    check.set_defaults(handler=check_command)

    devices = subparsers.add_parser("devices", help="list available audio devices")
    add_dry_run(devices)
    devices.set_defaults(handler=devices_command)

    geometry = subparsers.add_parser("geometry", help="generate a microphone coordinate file")
    geometry.add_argument("--layout", choices=("photographed", "custom"), default="photographed")
    geometry.add_argument("--spread", type=float, default=0.6, help="custom left/right UMA-16 distance in metres")
    geometry.add_argument("--back", type=float, default=0.3, help="custom UMA-16 distance behind UMA-8 in metres")
    geometry.add_argument("--height16", type=float, default=0.0, help="custom UMA-16 height relative to UMA-8 in metres")
    geometry.add_argument("--facing", choices=("up", "front"), default="front")
    geometry.add_argument("--out", type=Path, default=DEFAULT_POSITIONS)
    add_dry_run(geometry)
    geometry.set_defaults(handler=geometry_command)

    run = subparsers.add_parser("run", help="start unified real-time 2D direction estimation")
    run.add_argument("--devices", type=int, nargs=3, metavar=("A0", "A1", "A2"))
    run.add_argument("--channels", type=int, nargs=3, default=(8, 16, 16), metavar=("N0", "N1", "N2"))
    run.add_argument("--positions", type=Path)
    run.add_argument("--sync-channels", type=int, nargs=3, default=(0, 0, 0))
    run.add_argument("--speaker-position", type=float, nargs=2, default=(0.0, 1.0), metavar=("X", "Y"))
    run.add_argument("--max-sync-error", type=float)
    run.add_argument("--fs", type=int, default=48000)
    run.add_argument("--calibration-seconds", type=float, default=45.0)
    run.add_argument("--frame-ms", type=int, default=200)
    run.add_argument("--update-ms", type=int, default=500)
    run.add_argument("--coarse-step", type=float, default=0.10)
    run.add_argument("--refinements", type=int, default=3)
    run.add_argument("--min-score", type=float, default=0.25)
    run.add_argument(
        "--bounds", type=float, nargs=6, default=(-1, 1, -1, 1, 0, 0),
        metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"),
    )
    run.add_argument("--plot", action="store_true", help="show the live position plot")
    add_dry_run(run)
    run.set_defaults(handler=run_command)
    return parser


def interactive_main(parser: argparse.ArgumentParser) -> int:
    actions = {
        "1": ["check"],
        "2": ["devices"],
        "3": ["geometry"],
        "4": ["run", "--plot"],
    }
    while True:
        print("\n=== Triangulation Launcher ===")
        print("1. Check setup")
        print("2. List audio devices")
        print("3. Generate photographed microphone positions")
        print("4. Start unified 2D direction estimation")
        print("0. Exit")
        try:
            choice = input("Choose an option: ").strip()
            if choice == "0":
                return 0
            if choice not in actions:
                print("Choose a number from 0 to 4.")
                continue
            args = parser.parse_args(actions[choice])
            args.handler(args)
            input("\nPress Enter to return to the menu...")
        except (ValueError, OSError) as error:
            print(f"Error: {error}")
        except EOFError:
            print()
            return 0


def main() -> int:
    parser = build_parser()
    if len(sys.argv) == 1:
        return interactive_main(parser)
    args = parser.parse_args()
    if not hasattr(args, "handler"):
        parser.print_help()
        return 0
    try:
        return args.handler(args)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    except KeyboardInterrupt:
        print("\nStopped.")
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

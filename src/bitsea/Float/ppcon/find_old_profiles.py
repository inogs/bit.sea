#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


PROFILE_FILENAME = re.compile(r"SR(\d{7})_(\d{3})\.nc")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Report PPCON SR profiles that have a matching new SD profile (same WMO number and rise up number)."
    )
    parser.add_argument(
        "--inputdir",
        "-i",
        type=Path,
        required=True,
        help="Path to the directory: ONLINE/PPCON, ONLINE/CORIOLIS, ONLINE/SUPERFLOAT.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Remove matching SR profiles instead of reporting them.",
    )
    args = parser.parse_args()

    if not args.inputdir.is_dir():
        parser.error(f"not a directory: {args.inputdir}")

    for wmo_dir in args.inputdir.iterdir():
        if not wmo_dir.is_dir() or not re.fullmatch(r"\d{7}", wmo_dir.name):
            continue

        profiles_to_report = []
        for profile in wmo_dir.iterdir():
            match = PROFILE_FILENAME.fullmatch(profile.name)
            if not profile.is_file() or match is None:
                continue

            wmo_number, rise_number = match.groups()
            if wmo_number != wmo_dir.name:
                continue

            new_profile = wmo_dir / f"SD{wmo_number}_{rise_number}.nc"
            if new_profile.is_file():
                profiles_to_report.append((profile, wmo_number, rise_number))

        for profile, wmo_number, rise_number in profiles_to_report:
            if args.force:
                profile.unlink()
            else:
                print(f"SD/SR profiles pair exist for float {wmo_number} rise {rise_number}")


if __name__ == "__main__":
    main()

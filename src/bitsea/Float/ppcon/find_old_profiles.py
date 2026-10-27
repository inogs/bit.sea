#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


PROFILE_FILENAME = re.compile(r"SD(\d{7})_(\d{3})\.nc")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Report PPCON SD profiles that have a matching old SR profile."
    )
    parser.add_argument(
        "parent_dir",
        type=Path,
        help="Path to the ONLINE/PPCON directory.",
    )
    args = parser.parse_args()

    if not args.parent_dir.is_dir():
        parser.error(f"not a directory: {args.parent_dir}")

    for wmo_dir in args.parent_dir.iterdir():
        if not wmo_dir.is_dir() or not re.fullmatch(r"\d{7}", wmo_dir.name):
            continue

        for profile in wmo_dir.iterdir():
            match = PROFILE_FILENAME.fullmatch(profile.name)
            if not profile.is_file() or match is None:
                continue

            wmo_number, rise_number = match.groups()
            if wmo_number != wmo_dir.name:
                continue

            old_profile = wmo_dir / f"SR{wmo_number}_{rise_number}.nc"
            if old_profile.is_file():
                print(
                    f"old profile exists for: float n {wmo_number}, "
                    f"rise n {rise_number}"
                )


if __name__ == "__main__":
    main()

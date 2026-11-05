#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


PROFILE_FILENAME = re.compile(r"SR(\d{7})_(\d{3})\.nc")


def iter_wmo_directory_contents(inpdircontent: Path):
    root_dirs = []
    directory_contents = {}
    current_dir = None
    in_root_list = False

    for line in inpdircontent.read_text(encoding="utf-8").splitlines():
        if not line:
            continue

        if line == ".:":
            in_root_list = True
            current_dir = None
            continue

        if line.startswith("./") and line.endswith(":"):
            in_root_list = False
            current_dir = line[2:-1]
            if re.fullmatch(r"\d{7}", current_dir):
                directory_contents.setdefault(current_dir, [])
            else:
                current_dir = None
            continue

        if in_root_list:
            if re.fullmatch(r"\d{7}", line):
                root_dirs.append(line)
            continue

        if current_dir is not None:
            directory_contents.setdefault(current_dir, []).append(line)

    for wmo_dir in root_dirs:
        if wmo_dir in directory_contents:
            yield wmo_dir, directory_contents[wmo_dir]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Report PPCON SR profiles that have a matching new SD profile (same WMO number and rise up number)."
    )
    parser.add_argument(
        "--inputdir",
        "-I",
        type=Path,
        required=True,
        help="Path to the directory: ONLINE/PPCON, ONLINE/CORIOLIS, ONLINE/SUPERFLOAT.",
    )
    parser.add_argument(
        "--inpdircontent",
        "-i",
        type=Path,
        required=True,
        help="Path to the ls -1R text output for ONLINE/PPCON, ONLINE/CORIOLIS or ONLINE/SUPERFLOAT.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Remove matching SR profiles instead of reporting them.",
    )
    args = parser.parse_args()

    if not args.inputdir.is_dir():
        parser.error(f"not a directory: {args.inputdir}")

    if not args.inpdircontent.is_file():
        parser.error(f"not a file: {args.inpdircontent}")

    for wmo_dir, directory_contents in iter_wmo_directory_contents(args.inpdircontent):
        profiles_to_report = []
        for profile_name in directory_contents:
            if not profile_name.startswith("SR"):
                continue

            match = PROFILE_FILENAME.fullmatch(profile_name)
            if match is None:
                continue

            wmo_number, rise_number = match.groups()
            if wmo_number != wmo_dir:
                continue

            new_profile = f"SD{wmo_number}_{rise_number}.nc"
            if new_profile in directory_contents:
                profiles_to_report.append((profile_name, wmo_number, rise_number))

        for profile_name, wmo_number, rise_number in profiles_to_report:
            print(f"SD/SR profiles pair exist for float {wmo_number} rise {rise_number}")
            if args.force:
                removal_path = args.inputdir / wmo_number / profile_name
                print(f"Removing {removal_path}")
                if removal_path.is_file():
                    removal_path.unlink()

if __name__ == "__main__":
    main()

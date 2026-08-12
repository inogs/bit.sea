#!/usr/bin/env python3
import argparse
import os
import shutil
from pathlib import Path


def argument() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""Filter UPDATE_FILE_PATH rows and maintain PPCON copies from SUPERFLOAT.""")
    parser.add_argument(
        "--inputfile",
        "-i",
        required=True,
        help="Path to the update file (input list of float files).",
    )
    parser.add_argument(
        "--outputfile",
        "-o",
        required=True,
        help="Path to the local filtered file to write.",
    )
    parser.add_argument(
        "--inputdir",
        "-I",
        required=True,
        help="Path to the directory containing the input superfloat files.",
    )
    parser.add_argument(
        "--outputdir",
        "-O",
        required=True,
        help="Path to the directory where the PPCON files will be maintained.",
    )
    return parser.parse_args()


def _strip_prefix_once(value: str, prefix: str) -> str:
    return value.replace(prefix, "", 1)


def main():
    args = argument()

    update_file_path = Path(args.inputfile)
    local_input = Path(args.outputfile)
    superfloat_dir = Path(args.inputdir)
    ppcon_dir = Path(args.outputdir)

    local_input.parent.mkdir(parents=True, exist_ok=True)

    with update_file_path.open("r", encoding="utf-8") as src, local_input.open("w", encoding="utf-8") as dst:
        for raw_line in src:
            line = raw_line.rstrip("\n")
            if not line:
                continue

            filename = line.split(",", 1)[0]
            v1 = _strip_prefix_once(filename, "coriolis/")
            v2 = _strip_prefix_once(v1, "profiles/")
            superfloat_file = superfloat_dir / v2

            # Keep only lines that already exist in SUPERFLOAT.
            if superfloat_file.is_file():
                dst.write(line + "\n")

                # Maintain PPCON copy of the file.

                wmo = v2.split("/", 1)[0]
                ppcon_wmo = ppcon_dir / wmo
                ppcon_wmo.mkdir(parents=True, exist_ok=True)

                dst_file = ppcon_dir / v2
                dst_file.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(superfloat_file, dst_file)

if __name__ == "__main__":
   main()
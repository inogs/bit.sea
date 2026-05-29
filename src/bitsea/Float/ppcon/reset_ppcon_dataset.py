#!/usr/bin/env python3
import argparse
import os
import shutil
from pathlib import Path


def argument() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""Filter UPDATE_FILE_PATH rows and maintain PPCON copies from SUPERFLOAT.""")
    parser.add_argument(
        "--UPDATE_FILE_PATH",
        required=True,
        help="Path to the update file (input CSV-like text).",
    )
    parser.add_argument(
        "--LOCAL_INPUT",
        required=True,
        help="Path to the local filtered file to write.",
    )
    return parser.parse_args()


def _strip_prefix_once(value: str, prefix: str) -> str:
    return value.replace(prefix, "", 1)


def main() -> str:
    args = argument()

    online_repo = os.environ.get("ONLINE_REPO")
    if not online_repo:
        raise RuntimeError("ONLINE_REPO environment variable is not set.")

    update_file_path = Path(args.UPDATE_FILE_PATH)
    local_input = Path(args.LOCAL_INPUT)

    local_input.parent.mkdir(parents=True, exist_ok=True)

    with update_file_path.open("r", encoding="utf-8") as src, local_input.open("w", encoding="utf-8") as dst:
        for raw_line in src:
            line = raw_line.rstrip("\n")
            if not line:
                continue

            filename = line.split(",", 1)[0]
            v1 = _strip_prefix_once(filename, "coriolis/")
            v2 = _strip_prefix_once(v1, "profiles/")
            superfloat_file = Path(online_repo) / "SUPERFLOAT" / v2

            # Keep only lines that already exist in SUPERFLOAT.
            if superfloat_file.is_file():
                dst.write(line + "\n")

                wmo = v2.split("/", 1)[0]
                ppcon_wmo = Path(online_repo) / "PPCON" / wmo
                ppcon_wmo.mkdir(parents=True, exist_ok=True)

                dst_file = Path(online_repo) / "PPCON" / v2
                dst_file.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(superfloat_file, dst_file)

        nrows_int = sum(1 for _ in dst)
        nrows = str(nrows_int)

    return nrows


if __name__ == "__main__":
   nr_value = main()
   # Print only NR so shell can capture it with command substitution.
   print(nr_value)
import argparse
from pathlib import Path

from bitsea.commons.utils import file2stringlist


def argument():
    parser = argparse.ArgumentParser(description = '''
    Find differences between two argo floats files
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--input_NEW',"-N",
                                type = str,
                                required = True,
                                help = 'input file argo NEW, eg=')
    parser.add_argument(   '--input_OLD',"-O",
                                type = str,
                                required = True,
                                help = 'input file argo txt OLD')
    parser.add_argument(   '--outputfile',"-o",
                                type = str,
                                required = True,
                                help = 'output file name')

    return parser.parse_args()


def build_difference_file(input_new: Path | str, input_old: Path | str, output_file: Path | str) -> None:
    old_list = file2stringlist(str(input_old))
    new_list = file2stringlist(str(input_new))

    diff = []
    for line in new_list:
        if line not in old_list:
            diff.append(line)

    lines = []
    for line in diff:
        lines.append(line + "\n")

    with open(output_file, "wt") as file_handle:
        file_handle.writelines(lines)


def main() -> None:
    args = argument()
    build_difference_file(input_new=args.input_NEW, input_old=args.input_OLD, output_file=args.outputfile)


if __name__ == "__main__":
    main()


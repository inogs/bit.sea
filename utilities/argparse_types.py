"""
This module introduces some function that can be used as "types" arguments
while configuring argparse.
"""
import argparse
from datetime import datetime
from pathlib import Path


def generic_path(arg_path: str) -> Path:
    return Path(arg_path)


def existing_dir_path(arg_path: str) -> Path:
    arg_path = generic_path(arg_path)

    if not arg_path.exists():
        raise argparse.ArgumentTypeError(
            'Path "{}" does not exist.'.format(arg_path)
        )

    if not arg_path.is_dir():
        raise argparse.ArgumentTypeError(
            'Path "{}" is not a directory.'.format(arg_path)
        )

    return arg_path


def existing_file_path(arg_path: str) -> Path:
    arg_path = generic_path(arg_path)

    if not arg_path.exists():
        raise argparse.ArgumentTypeError(
            'Path "{}" does not exist.'.format(arg_path)
        )

    if not arg_path.is_file():
        raise argparse.ArgumentTypeError(
            'Path "{}" does not point to a file.'.format(arg_path)
        )

    return arg_path


def path_inside_an_existing_dir(arg_path: str) -> Path:
    arg_path = generic_path(arg_path).resolve()

    if not arg_path.parent.exists():
        raise argparse.ArgumentTypeError(
            'Path "{}" does not exist.'.format(arg_path)
        )

    return arg_path


def non_existing_path_inside_an_existing_dir(arg_path: str) -> Path:
    arg_path = path_inside_an_existing_dir(arg_path)

    if arg_path.exists():
        raise argparse.ArgumentTypeError(
            'Path "{}" already exists'.format(arg_path)
        )

    return arg_path


def date_from_str(arg_path: str) -> datetime:
    formats = ('%Y-%m-%d', '%Y-%m-%dT%H:%M', '%Y-%m-%dT%H:%M:%S')

    output_date = None
    for date_format in formats:
        try:
            output_date = datetime.strptime(arg_path, date_format)
        except ValueError:
            continue
        break

    if output_date is None:
        raise argparse.ArgumentTypeError(
            'Argument "{}" can not be interpreted as a date; valid formats '
            'are: {}'.format(
                arg_path,
                ','.join(('"{}"'.format(f) for f in formats))
            )
        )

    return output_date

def some_among(choices):
    """
    Arguments:
    choices : list of strings
    Returns: function

    """
    def check_if_valids(inputstr):
        strlist=inputstr.split(",")
        for entry_name in strlist:
            if entry_name not in choices:
                raise argparse.ArgumentTypeError(entry_name + " is not a valid choice. Valid choice are " + ", ".join(choices))
        return strlist
    return check_if_valids

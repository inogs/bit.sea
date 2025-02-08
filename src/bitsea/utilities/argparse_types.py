"""
This module introduces some function that can be used as "types" arguments
while configuring argparse.
"""

import argparse
from collections import OrderedDict
from datetime import datetime
from datetime import timedelta
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


def dir_inside_an_existing_path(arg_path: str) -> Path:
    """
    This function should be used as a type in argparse when the argument being
    parsed might be an existing directory or a directory that has yet to be created.
    If the directory does not already exist, this function will create it before
    returning a Path object that refers to it.

    Parameters:
        arg_path (str): The input path to be validated or created as a directory
        if it does not exist.

    Returns:
        Path: A Path object that refers to the validated or created directory.

    Raises:
        argparse.ArgumentTypeError: If the argument refers to a directory that
        cannot be created due to non-existent parent directory.

        argparse.ArgumentTypeError: If the argument refers to a directory that
        cannot be created because its proposed parent is not a directory.

        argparse.ArgumentTypeError: If the argument refers to a location that
        exists but is not a directory.

    """

    arg_path = generic_path(arg_path)
    if not arg_path.parent.exists():
        raise argparse.ArgumentTypeError(
            'Path "{}" can not be created because "{}" does not exist.'.format(
                arg_path, arg_path.parent
            )
        )
    if not arg_path.parent.is_dir():
        raise argparse.ArgumentTypeError(
            'Path "{}" can not be created because "{}" is not a '
            "directory.".format(arg_path, arg_path.parent)
        )
    if not arg_path.exists():
        arg_path.mkdir()
    if not arg_path.is_dir():
        raise argparse.ArgumentTypeError(
            'Path "{}" is not a directory.'.format(arg_path)
        )
    return arg_path


def amount_of_time(default_unit="s"):
    # I used an OrderedDict only to be sure that when I print the allowed
    # chars in the error message (or in the help) I see them in the right
    # order
    allowed_units = OrderedDict(
        (
            ("s", 1),
            ("m", 60),
            ("h", 60 * 60),
            ("d", 60 * 60 * 24),
        )
    )

    def read_amount_of_time(t: str):
        with_suffix = False
        try:
            t = int(t)
        except ValueError:
            with_suffix = True

        if with_suffix:
            t_value_str = t[:-1]
            t_suffix = t[-1]
        else:
            t_value_str = t
            t_suffix = default_unit

        invalid_entry = False
        if t_suffix not in allowed_units:
            invalid_entry = True

        t_value = 0
        try:
            t_value = int(t_value_str)
        except ValueError:
            invalid_entry = True

        if invalid_entry:
            raise argparse.ArgumentTypeError(
                "Invalid amount of time specified. It must be an integer "
                "number followed by one of the following letters: {}. "
                'Received: "{}"'.format(
                    ", ".join(["{}".format(k) for k in allowed_units.keys()]), t
                )
            )

        return timedelta(seconds=t_value * allowed_units[t_suffix])

    return read_amount_of_time


def date_from_str(arg_path: str) -> datetime:
    formats = (
        "%Y%m%d",
        "%Y%m%d-%H:%M",
        "%Y%m%d-%H:%M:%S",
        "%Y-%m-%d",
        "%Y-%m-%dT%H:%M",
        "%Y-%m-%dT%H:%M:%S",
    )

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
            "are: {}".format(
                arg_path, ",".join(('"{}"'.format(f) for f in formats))
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
        strlist = inputstr.split(",")
        for entry_name in strlist:
            if entry_name not in choices:
                raise argparse.ArgumentTypeError(
                    entry_name
                    + " is not a valid choice. Valid choice are "
                    + ", ".join(choices)
                )
        return strlist

    return check_if_valids

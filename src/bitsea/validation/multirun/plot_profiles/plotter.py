from argparse import ArgumentParser
from pathlib import Path
from warnings import warn

from bitsea.basins import V2
from bitsea.utilities.argparse_types import existing_file_path
from bitsea.validation.multirun.plot_profiles import draw_profile_plots
from bitsea.validation.multirun.plot_profiles.tools.read_config import (
    read_config_from_file,
)
from bitsea.validation.multirun.plot_profiles.tools.read_config import (
    read_output_dir,
)


try:
    import mpi4py.MPI
except Exception as e:
    warn(
        "mpi4py can not be imported. This code will run serially. The reason "
        f"of the error was:\n{e}"
    )


MAIN_DIR = Path(__file__).resolve().parent
CONFIG_FILE = MAIN_DIR / "config.yaml"


BASINS = tuple(V2.P.basin_list)


def configure_argparse():
    parser = ArgumentParser()
    parser.description = (
        "The MultirunProfilePlotter is a script to plot temporal series and "
        "depth profiles of several different runs of one or more models. "
        "You may configure the behaviour of the script by setting the values "
        "of the config.yaml file in the main directory of this script or "
        "writing another configuration file yourself."
    )
    parser.add_argument(
        "config_file",
        type=existing_file_path,
        nargs="?",
        default=CONFIG_FILE,
        help="The path of the config file used by this script. By default, it "
        "uses {}".format(CONFIG_FILE),
    )

    parser.add_argument(
        "--output_dir",
        "-o",
        type=read_output_dir,
        default=None,
        help="The path where this script will save the output plots",
    )
    return parser.parse_args()


def main():
    args = configure_argparse()
    config = read_config_from_file(args.config_file, output_dir=args.output_dir)

    draw_profile_plots(config, BASINS)


if __name__ == "__main__":
    main()

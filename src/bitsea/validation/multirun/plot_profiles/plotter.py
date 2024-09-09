from argparse import ArgumentParser
from pathlib import Path

from bitsea.basins import V2

from bitsea.validation.multirun.plot_profiles.tools.read_config import \
    InvalidConfigFile, read_config_from_file, read_output_dir
from bitsea.validation.multirun.plot_profiles import draw_profile_plots


try:
    import mpi4py
    isParallel = True
except ModuleNotFoundError:
    isParallel = False


MAIN_DIR = Path(__file__).resolve().parent
CONFIG_FILE = MAIN_DIR / 'config.yaml'


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
        'config_file',
        type=str,
        nargs='?',
        default=CONFIG_FILE,
        help='The path of the config file used by this script. By default, it '
             'uses {}'.format(CONFIG_FILE)
    )

    parser.add_argument(
        '--output_dir',
        '-o',
        type=read_output_dir,
        default=None,
        help='The path where this script will save the output plots'
    )
    return parser.parse_args()


def main():
    args = configure_argparse()
    config = read_config_from_file(args.config_file)

    if args.output_dir is None and config.output_options.output_dir is None:
        raise InvalidConfigFile(
            'output_dir has not been specified in the output section of the '
            'config file. Specify an output directory in the config file or '
            'use the --output-dir flag from the command line.'
        )

    draw_profile_plots(config, BASINS, args.output_dir)


if __name__ == '__main__':
    main()

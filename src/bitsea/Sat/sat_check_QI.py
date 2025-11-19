import argparse
import logging
from collections.abc import Sequence
from contextlib import ExitStack
from datetime import datetime
from pathlib import Path
from sys import exit as sys_exit
from typing import List

import netCDF4
import numpy as np

from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.postproc import masks
from bitsea.Sat import SatManager as Sat
from bitsea.utilities.argparse_types import existing_dir_path
from bitsea.utilities.argparse_types import some_among
from bitsea.utilities.mpi_serial_interface import get_mpi_communicator


if __name__ == "__main__":
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser(
        description="""
            Apply check based on QI provided by OCTAC
            Produces CHECKED files for each date at satellite resolution.
            QI = (value - ClimatologyMedianData)/ClimatologyStandardDeviation

            """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--origdir",
        "-i",
        type=existing_dir_path,
        required=True,
        help="ORIG sat directory",
    )

    parser.add_argument(
        "--checkdir",
        "-o",
        type=existing_dir_path,
        required=True,
        help="Base for CHECKED sat directory",
    )

    parser.add_argument(
        "--mesh",
        "-m",
        type=str,
        required=True,
        choices=[
            "SatOrigMesh",
            "V4mesh",
            "V1mesh",
            "KD490mesh",
            "SAT1km_mesh",
            "Mesh24",
        ],
        help="Name of the mesh of sat ORIG and used to dump checked data.",
    )
    parser.add_argument(
        "--varnames",
        "-v",
        type=some_among(
            [
                "CHL",
                "KD490",
                "DIATO",
                "NANO",
                "PICO",
                "DINO",
                "RRS412",
                "RRS443",
                "RRS490",
                "RRS510",
                "RRS555",
                "RRS670",
            ]
        ),
        required=True,
        help="Var name, corresponding to P_l, kd490, P1l, P2l P3l, P4l",
    )
    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Overwrite existing variables in files",
    )
    parser.add_argument(
        "--QI",
        required=True,
        type=str,
        help="QI threshold, usually 2, 2.5 or 3",
    )
    parser.add_argument(
        "--Kd_min",
        required=False,
        type=str,
        default="0.021",
        help="Kd minimum value threshold",
    )
    parser.add_argument("--serial", action="store_true", help="Do not use mpi")

    return parser.parse_args()


def sat_check(
    *,
    inputfiles: Sequence[Path],
    outputdir: Path,
    mesh: str,
    varnames: List[str],
    qi_threshold: float,
    Kd_min: float = 0.021,
    force: bool = False,
    mpi_communicator=None,
) -> List[Path]:
    if mpi_communicator is None:
        mpi_communicator = get_mpi_communicator()

    rank = mpi_communicator.Get_rank()
    nranks = mpi_communicator.size

    maskSat = getattr(masks, mesh)

    checked_file_list = [outputdir / f.name for f in inputfiles]

    owned_files = inputfiles[rank::nranks]
    LOGGER.info("Process rank %i will update %i files", rank, len(owned_files))

    for filename in owned_files:
        outfile = outputdir / filename.name
        LOGGER.debug("Checking file %s...", outfile)

        input_file = None
        with ExitStack() as data_files:
            for varname in varnames:
                writing_mode = "a" if outfile.exists() else "w"

                output_file = data_files.enter_context(
                    netCDF4.Dataset(outfile, writing_mode)
                )

                if writing_mode == "w":
                    Sat.write_coords_in_sat_file(output_file, mesh=maskSat)

                available_variables = set(output_file.variables.keys())

                condition_to_write = varname not in available_variables
                if force:
                    condition_to_write = True

                if not condition_to_write:
                    LOGGER.debug(
                        "File %s already contains variable %s", outfile, varname
                    )
                    continue

                LOGGER.debug(
                    "We need to process variable %s; reading its content "
                    "from %s",
                    varname,
                    filename,
                )
                if input_file is None:
                    LOGGER.debug("Opening file %s", filename)
                    input_file = data_files.enter_context(
                        netCDF4.Dataset(filename, "r")
                    )
                else:
                    LOGGER.debug("File %s was already open", filename)

                qi_var_name = f"QI_{varname}"
                if varname in ["DIATO", "NANO", "PICO", "DINO"]:
                    qi_var_name = "QI_CHL"

                LOGGER.debug(
                    "Reading variable %s from %s", qi_var_name, filename
                )
                sat_qi = np.asarray(input_file.variables[qi_var_name][:])
                LOGGER.debug("Reading variable %s from %s", varname, filename)
                sat_values = np.asarray(input_file.variables[varname][:])

                if len(sat_qi.shape) == 3:
                    sat_qi = sat_qi[0, :, :]
                if len(sat_values.shape) == 3:
                    sat_values = sat_values[0, :, :]

                if varname == "KD490":
                    # sat_values[(sat_values<Kd_min) & (sat_values>0)] = Kd_min
                    sat_values[(sat_values < Kd_min) & (sat_values > 0)] = (
                        Sat.fillValue
                    )  # Filter out values below Kd490_min threshold

                bad = np.abs(sat_qi) > qi_threshold  # 2.0
                sat_values[bad] = Sat.fillValue

                LOGGER.info(
                    "Writing variable %s inside file %s", outfile, varname
                )
                if varname not in available_variables:
                    var_data = output_file.createVariable(
                        varname,
                        "f",
                        ("time", "lat", "lon"),
                        zlib=True,
                        shuffle=True,
                        complevel=4,
                        fill_value=Sat.fillValue,
                    )
                    setattr(var_data, "missing_value", Sat.fillValue)
                    setattr(var_data, "fillValue", Sat.fillValue)
                else:
                    var_data = output_file.variables[varname]

                d = datetime.now()
                setattr(
                    var_data, "creation_time", d.strftime("%Y%m%d-%H:%M:%S")
                )
                var_data[:] = sat_values

    return checked_file_list


def configure_logger(rank: int = 0) -> None:
    format_string = (
        f"%(asctime)s [rank={rank:0>3}] - %(name)s - "
        "%(levelname)s - %(message)s"
    )
    formatter = logging.Formatter(format_string)

    LOGGER.setLevel(logging.INFO)

    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(formatter)

    LOGGER.addHandler(handler)


def main():
    args = argument()

    if not args.serial:
        # noinspection PyUnresolvedReferences
        import mpi4py.MPI

    mpi_communicator = get_mpi_communicator()
    configure_logger(rank=mpi_communicator.Get_rank())

    time_start = "19501231"
    time_end = "20500101"
    TI = TimeInterval(time_start, time_end, "%Y%m%d")
    TL_orig = TimeList.fromfilenames(
        TI, args.origdir, "*.nc", prefix="", dateformat="%Y%m%d"
    )

    sat_check(
        inputfiles=TL_orig.filelist,
        outputdir=args.checkdir,
        mesh=args.mesh,
        varnames=args.varnames,
        qi_threshold=float(args.QI),
        Kd_min=float(args.Kd_min),
        force=args.force,
        mpi_communicator=mpi_communicator,
    )
    return 0


if __name__ == "__main__":
    sys_exit(main())

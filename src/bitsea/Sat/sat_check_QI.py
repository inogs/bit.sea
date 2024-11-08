import argparse
from pathlib import Path
from sys import exit as sys_exit
from typing import Iterable
from typing import List

import numpy as np
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.postproc import masks
from bitsea.Sat import SatManager as Sat
from bitsea.utilities.argparse_types import existing_dir_path
from bitsea.utilities.argparse_types import some_among
from bitsea.utilities.mpi_serial_interface import get_mpi_communicator


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


args = argument()


def sat_check(
    *,
    inputfiles: Iterable[Path],
    outputdir: Path,
    mesh: str,
    varnames: List[str],
    QI: float,
    Kd_min: float = float(0.021),
    force: bool = False,
    serial: bool = False,
) -> List[Path]:
    if not serial:
        pass

    comm = get_mpi_communicator()
    rank = comm.Get_rank()
    nranks = comm.size

    CHECKDIR = outputdir
    THRESHOLD = QI

    maskSat = getattr(masks, mesh)

    checkedFileList = [CHECKDIR / f.name for f in inputfiles]

    for filename in inputfiles[rank::nranks]:
        outfile = CHECKDIR / filename.name

        for varname in varnames:
            writing_mode = Sat.writing_mode(outfile)

            condition_to_write = not Sat.exist_valid_variable(varname, outfile)
            if force:
                condition_to_write = True
            if not condition_to_write:
                continue

            if varname in ["DIATO", "NANO", "PICO", "DINO"]:
                QI = Sat.readfromfile(filename, "QI_CHL")
            else:
                QI = Sat.readfromfile(filename, "QI_" + varname)
            VALUES = Sat.readfromfile(filename, varname)

            if varname == "KD490":
                # VALUES[(VALUES<Kd_min) & (VALUES>0)] = Kd_min
                VALUES[(VALUES < Kd_min) & (VALUES > 0)] = (
                    Sat.fillValue
                )  # Filter out values below Kd490_min threshold

            bad = np.abs(QI) > THRESHOLD  # 2.0
            VALUES[bad] = Sat.fillValue
            print(outfile, varname, flush=True)
            Sat.dumpGenericfile(
                outfile, VALUES, varname, mesh=maskSat, mode=writing_mode
            )
    return checkedFileList


def main():
    args = argument()

    Timestart = "19501231"
    Time__end = "20500101"
    TI = TimeInterval(Timestart, Time__end, "%Y%m%d")
    TL_orig = TimeList.fromfilenames(
        TI, args.origdir, "*.nc", prefix="", dateformat="%Y%m%d"
    )

    sat_check(
        inputfiles=TL_orig.filelist,
        outputdir=args.checkdir,
        mesh=args.mesh,
        varnames=args.varnames,
        QI=float(args.QI),
        Kd_min=float(args.Kd_min),
        force=args.force,
        serial=args.serial,
    )
    return 0


if __name__ == "__main__":
    sys_exit(main())

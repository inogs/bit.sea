import argparse
import logging
from datetime import datetime
from pathlib import Path
from sys import exit as sys_exit
from typing import List
from typing import Sequence

import bitsea.Sat.SatManager as Sat
import numpy as np
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.postproc import masks
from bitsea.utilities.argparse_types import date_from_str
from bitsea.utilities.argparse_types import existing_dir_path
from bitsea.utilities.argparse_types import some_among
from bitsea.utilities.mpi_serial_interface import get_mpi_communicator


def configure_argparse():
    parser = argparse.ArgumentParser(
        description="""
    Generic averager for sat files.
    It works with one mesh, without performing interpolations.
    Empty dirdates is provided when the number of input daily files is lesser than 3.
    """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--checkdir",
        "-i",
        type=existing_dir_path,
        required=True,
        help=""" CHECKED sat directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/""",
    )
    parser.add_argument(
        "--outdir",
        "-o",
        type=existing_dir_path,
        required=True,
        help=""" OUT average dates sat directory""",
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
            "Mesh4",
        ],
        help=""" Name of the mesh of sat ORIG and used to dump checked data.""",
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
    )

    parser.add_argument(
        "--timeaverage",
        "-t",
        type=str,
        required=True,
        choices=[
            "monthly",
            "weekly_tuesday",
            "weekly_friday",
            "weekly_monday",
            "weekly_thursday",
            "tendays",
        ],
        help=""" Name of the mesh of sat ORIG and used to dump checked data.""",
    )
    parser.add_argument(
        "--ignore-after",
        "-a",
        type=date_from_str,
        required=False,
        default=None,
        help="Ignore all input files with dates later than the one submitted here",
    )
    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="""Overwrite existing variables in files
                                """,
    )
    parser.add_argument(
        "--serial", action="store_true", help="""Do not use mpi"""
    )
    return parser.parse_args()


def aveSat(
    *,
    inputfiles: Sequence[Path],
    datetimes: List[datetime],
    outputdir: Path,
    mesh: str,
    varnames: List[str],
    outfrequency: str = "weekly",
    week_day: int = 1,
    force: bool = False,
) -> List[Path]:
    valid_frequencies = ("monthly", "weekly", "ten_days")
    if outfrequency not in valid_frequencies:
        raise ValueError(
            f"outfrequency must be one among {valid_frequencies}; received {outfrequency}"
        )
    if week_day < 0 or week_day >= 7:
        raise ValueError(
            "week_day must be greater or equal than 0 and "
            f"smaller than 7; received {week_day}"
        )
    if int(week_day) != week_day:
        raise ValueError(
            f"week_day must be an integer number: received {week_day}"
        )

    comm = get_mpi_communicator()
    rank = comm.Get_rank()
    nranks = comm.size

    logger = logging.getLogger(__name__)

    mask_sat = getattr(masks, mesh)

    suffix = inputfiles[0].name[8:]
    TLCheck = TimeList(datetimes)
    if outfrequency == "monthly":
        TIME_reqs = TLCheck.getMonthlist()
    elif outfrequency == "weekly":
        TIME_reqs = TLCheck.getWeeklyList(week_day)
    elif outfrequency == "ten_days":
        TIME_reqs = TLCheck.getSpecificIntervalList(10, "19971001-12:00:00")
    else:
        raise ValueError(
            f'outfrequency must be one among "monthly", "weekly" or "ten_days"; received {outfrequency}'
        )

    jpi = mask_sat.jpi
    jpj = mask_sat.jpj
    averaged_filelist = [outputdir / (req.string + suffix) for req in TIME_reqs]

    for req in TIME_reqs[rank::nranks]:
        outfile = outputdir / (req.string + suffix)

        ii, w = TLCheck.select(req)
        n_files = len(ii)
        if n_files < 3:
            logger.info(
                "For %s there are less than 3 files - Skipping "
                "average generation; %s",
                req,
                " ".join([TLCheck.Timelist[k].strftime("%Y%m%d") for k in ii]),
            )
            continue

        for varname in varnames:
            writing_mode = Sat.writing_mode(outfile)
            if outfile.exists():
                try:
                    dateweek_string = Sat.read_variable_attribute(
                        outfile, varname, "average_of"
                    )
                except (AttributeError, IndexError):
                    dateweek_string = ""
                n_dates = len(dateweek_string.split(","))

                if n_files > n_dates:
                    logger.debug(
                        "For req %s there are more files (%s) than "
                        "dates (%s)",
                        req,
                        n_files,
                        n_dates,
                    )
                    pass
                else:
                    condition_to_write = not Sat.exist_valid_variable(
                        varname, outfile
                    )
                    if force:
                        condition_to_write = True
                    if not condition_to_write:
                        continue

            logger.info("Writing variable %s inside file %s", varname, outfile)
            dateweek = []

            M = np.zeros((n_files, jpj, jpi), np.float32)
            for iFrame, j in enumerate(ii):
                inputfile = inputfiles[j]
                frame_values = Sat.readfromfile(inputfile, varname)
                M[iFrame, :, :] = frame_values
                idate = datetimes[j]
                date8 = idate.strftime("%Y%m%d")
                dateweek.append(date8)
            if varname == "KD490":
                output_data = Sat.averager(M)
            else:
                output_data = Sat.logAverager(M)
            dateweek_string = ",".join(dateweek)
            var_attributes = {"average_of": dateweek_string}
            Sat.dumpGenericfile(
                outfile,
                output_data,
                varname,
                mesh=mask_sat,
                mode=writing_mode,
                var_attributes=var_attributes,
            )
    return averaged_filelist


def main():
    args = configure_argparse()

    if not args.serial:
        # noinspection PyUnresolvedReferences
        import mpi4py.MPI  # unused import

    logging.basicConfig(level=logging.INFO)

    time_start = datetime.strptime("19500101", "%Y%m%d")
    time_end = datetime.strptime("20500101", "%Y%m%d")

    if args.ignore_after is not None:
        time_end = args.ignore_after

    TI = TimeInterval.fromdatetimes(time_start, time_end)
    TL = TimeList.fromfilenames(
        TI, args.checkdir, "*.nc", prefix="", dateformat="%Y%m%d"
    )

    out_frequency = "weekly"
    weekday = None
    if args.timeaverage == "monthly":
        out_frequency = "monthly"
    if args.timeaverage == "weekly_tuesday":
        weekday = 2
    if args.timeaverage == "weekly_friday":
        weekday = 5
    if args.timeaverage == "weekly_monday":
        weekday = 1
    if args.timeaverage == "weekly_thursday":
        weekday = 4
    if args.timeaverage == "tendays":
        out_frequency = "tendays"

    aveSat(
        inputfiles=TL.filelist,
        datetimes=TL.Timelist,
        outputdir=args.outdir,
        mesh=args.mesh,
        varnames=args.varnames,
        outfrequency=out_frequency,
        week_day=weekday,
        force=args.force,
    )
    return 0


if __name__ == "__main__":
    sys_exit(main())

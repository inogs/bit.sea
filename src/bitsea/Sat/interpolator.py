import argparse

from bitsea.utilities.argparse_types import some_among


def argument():
    parser = argparse.ArgumentParser(
        description="""
    Interpolates from a fine mesh to a coarser output mesh.
    Works in parallel
    """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--inputdir",
        "-i",
        type=str,
        required=True,
        help=""" E.g. dir with files on 1km mesh (ORIG)""",
    )

    parser.add_argument(
        "--outputdir",
        "-o",
        type=str,
        required=True,
        help=""" E.g. dir with files on 1/24 mesh (interpolated)""",
    )
    parser.add_argument(
        "--inmesh",
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
        help=""" Name of the mesh of the input sat file""",
    )
    parser.add_argument(
        "--maskfile",
        "-m",
        type=str,
        required=True,
        help=""" Path of the meshmask corresponding to output sat files""",
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


args = argument()


from bitsea.commons.Timelist import TimeList
from bitsea.commons.time_interval import TimeInterval
from bitsea.Sat import SatManager as Sat
from bitsea.Sat import interp2d
from bitsea.commons.mask import Mask
from bitsea.postproc import masks
from bitsea.commons.utils import addsep
import os
from bitsea.utilities.mpi_serial_interface import get_mpi_communicator

from datetime import datetime

maskIn = getattr(masks, args.inmesh)


if not args.serial:
    pass


comm = get_mpi_communicator()
rank = comm.Get_rank()
nranks = comm.size

TheMask = Mask(args.maskfile)
x = TheMask.xlevels[0, :]
y = TheMask.ylevels[:, 0]

xOrig = maskIn.lon
yOrig = maskIn.lat

I_START, I_END = interp2d.array_of_indices_for_slicing(x, xOrig)
J_START, J_END = interp2d.array_of_indices_for_slicing(y, yOrig)

INPUTDIR = addsep(args.inputdir)
OUTPUTDIR = addsep(args.outputdir)
dateformat = "%Y%m%d"

Timestart = "19500101"
Time__end = "20500101"

TI = TimeInterval(Timestart, Time__end, "%Y%m%d")
TL = TimeList.fromfilenames(
    TI, INPUTDIR, "*.nc", prefix="", dateformat=dateformat
)


for filename in TL.filelist[rank::nranks]:
    outfile = OUTPUTDIR + os.path.basename(filename)
    condition_for_points = False

    for varname in args.varnames:
        writing_mode = Sat.writing_mode(outfile)

        attributes = {}
        try:
            attributes["average_of"] = Sat.read_variable_attribute(
                filename, varname, "average_of"
            )
        except AttributeError:
            pass
        try:
            attributes["from_variable_created_at"] = (
                Sat.read_variable_attribute(filename, varname, "creation_time")
            )
        except AttributeError:
            pass

        # condition_to write section ######
        condition_to_write = not Sat.exist_valid_variable(varname, outfile)

        if args.force:
            condition_to_write = True

        if not condition_to_write:
            try:
                attr_value = Sat.read_variable_attribute(
                    outfile, varname, "from_variable_created_at"
                )
            except AttributeError:
                attr_value = None

            if (
                attr_value is not None
            ) and "from_variable_created_at" in attributes:
                d_out = datetime.strptime(attr_value, "%Y%m%d-%H:%M:%S")
                d_in = datetime.strptime(
                    attributes["from_variable_created_at"], "%Y%m%d-%H:%M:%S"
                )
                condition_to_write = d_in > d_out
            else:
                condition_to_write = True

        if not condition_to_write:
            continue
            ##################################

        condition_for_points = True
        Mfine = Sat.readfromfile(filename, varname)
        print(outfile, varname, flush=True)

        Mout, usedPoints = interp2d.interp_2d_by_cells_slices(
            Mfine,
            TheMask,
            I_START,
            I_END,
            J_START,
            J_END,
            min_cov=0.0,
            ave_func=Sat.mean,
        )
        Sat.dumpGenericfile(
            outfile,
            Mout,
            varname,
            mesh=TheMask,
            mode=writing_mode,
            var_attributes=attributes,
        )

    if condition_for_points:
        Sat.dumpGenericfile(
            outfile, usedPoints, "Points", mesh=TheMask, mode="a"
        )

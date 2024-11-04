import argparse
from bitsea.utilities.argparse_types import some_among
from bitsea.utilities.argparse_types import existing_dir_path


def argument():
    parser = argparse.ArgumentParser(description = '''
    Apply check based on QI provided by OCTAC
    Produces CHECKED files for each date at satellite resolution.
    QI = (value - ClimatologyMedianData)/ClimatologyStandardDeviation

    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--origdir', '-i',
                                type = existing_dir_path,
                                required = True,
                                help = ''' ORIG sat directory'''
                                )

    parser.add_argument(   '--checkdir', '-o',
                                type = existing_dir_path,
                                required = True,
                                help = ''' Base for CHECKED sat directory'''
                                )

    parser.add_argument(   '--mesh', '-m',
                                type = str,
                                required = True,
                                choices = ['SatOrigMesh','V4mesh','V1mesh','KD490mesh','SAT1km_mesh', 'Mesh24'],
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
                                )
    parser.add_argument(   '--varnames', '-v',
                                type = some_among(['CHL','KD490','DIATO','NANO','PICO', 'DINO','RRS412','RRS443','RRS490','RRS510','RRS555','RRS670']),
                                required = True,
                                help = '''Var name, corresponding to P_l, kd490, P1l, P2l P3l, P4l'''
                                )
    parser.add_argument(   '--force', '-f',
                                action='store_true',
                                help = """Overwrite existing variables in files
                                """)
    parser.add_argument(   '--QI',
                                required=True,
                                type=str,
                                help = """QI threshold, usually 2, 2.5 or 3
                                """)
    parser.add_argument(   '--Kd_min',
                                required=False,
                                type=str,
                                default="0.021",
                                help = '''Kd minimum value threshold''')
    parser.add_argument(   '--serial',
                                action='store_true',
                                help = '''Do not use mpi''')


    return parser.parse_args()


args = argument()

from pathlib import Path
from typing import List
from bitsea.commons.Timelist import TimeList
from bitsea.commons.time_interval import TimeInterval
from bitsea.postproc import masks
import numpy as np
from bitsea.Sat import SatManager as Sat
from bitsea.utilities.mpi_serial_interface import get_mpi_communicator


def sat_check(*,
inputdir : Path,
outputdir: Path,
mesh: str,
varnames : List[str],
QI : float,
Kd_min: float = float(0.021),
force : bool = False,
serial : bool = False,
):
    if not serial:
        import mpi4py


    comm = get_mpi_communicator()
    rank  = comm.Get_rank()
    nranks =comm.size


    ORIGDIR = inputdir
    CHECKDIR = outputdir
    THRESHOLD = QI

    maskSat = getattr(masks,mesh)

    Timestart="19501231"
    Time__end="20500101"
    TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
    TL_orig = TimeList.fromfilenames(TI, ORIGDIR ,"*.nc",prefix='',dateformat='%Y%m%d')

    for filename in TL_orig.filelist[rank::nranks]:

        outfile = CHECKDIR / filename.name
        for varname in varnames:

            writing_mode = Sat.writing_mode(outfile)

            condition_to_write = not Sat.exist_valid_variable(varname,outfile)
            if args.force: condition_to_write=True
            if not condition_to_write: continue


            if varname in ["DIATO","NANO","PICO","DINO"]:
                QI = Sat.readfromfile(filename, "QI_CHL")
            else:
                QI = Sat.readfromfile(filename, "QI_" + varname)
            VALUES = Sat.readfromfile(filename,varname)

            if varname == 'KD490':
                # VALUES[(VALUES<Kd_min) & (VALUES>0)] = Kd_min
                VALUES[(VALUES<Kd_min) & (VALUES>0)] = Sat.fillValue   # Filter out values below Kd490_min threshold

            bad = np.abs(QI) > THRESHOLD # 2.0
            VALUES[bad] = Sat.fillValue
            print(outfile, varname, flush=True)
            Sat.dumpGenericfile(outfile, VALUES, varname, mesh=maskSat,mode=writing_mode)
    return 0

if __name__ == "__main__":
    args = argument()
    exit(
        sat_check(
            inputdir=args.origdir,
            outputdir=args.checkdir,
            mesh=args.mesh,
            varnames=args.varnames,
            QI = float(args.QI),
            Kd_min = float(args.Kd_min),
            force = args.force,
            serial=args.serial,
        )
    )


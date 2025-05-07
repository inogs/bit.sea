import argparse

import mpi4py.MPI
import numpy as np

import bitsea.matchup.matchup as matchup
import bitsea.Sat.SatManager as Satmodule
from bitsea.commons.dataextractor import DataExtractor
from bitsea.commons.layer import Layer
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.instruments.var_conversions import SAT_VARS
from bitsea.layer_integral.mapbuilder import MapBuilder
from bitsea.utilities.argparse_types import date_from_str
from bitsea.utilities.argparse_types import existing_dir_path
from bitsea.utilities.argparse_types import existing_file_path
from bitsea.utilities.argparse_types import some_among
from bitsea.utilities.mpi_serial_interface import get_mpi_communicator
from bitsea.validation.deliverables import netcdf_validation_file


def argument():
    parser = argparse.ArgumentParser(
        description="""
    Calculates Chl of Kd statistics on matchups of Model and Sat.
    Main result are
    - BGC_CLASS4_CHL_RMS_SURF_BASIN
    - BGC_CLASS4_CHL_BIAS_SURF_BASIN
    of deliverable CMEMS-Med-biogeochemistry-ScCP-1.0.pdf

    Other similar results are also calculated and dumped in
    netcdf files.
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--satdir",
        "-s",
        type=existing_dir_path,
        required=False,
        default="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI/MONTHLY_V4/",
        help=""" Input satellite dir""",
    )
    parser.add_argument(
        "--inputmodeldir",
        "-i",
        type=existing_dir_path,
        required=True,
        help=""" Input model dir, where P_l files are, usually ../wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/""",
    )
    parser.add_argument(
        "--outdir",
        "-o",
        type=existing_dir_path,
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--maskfile",
        "-m",
        type=existing_file_path,
        required=True,
        help="Path of the mask file",
    )
    parser.add_argument(
        "--coastness",
        "-c",
        type=some_among(["coast", "open_sea", "everywhere"]),
        required=True,
        help="definition of mask to apply to the statistics",
    )
    parser.add_argument(
        "--layer",
        "-l",
        type=int,
        required=True,
        help="Layer of the model, in meters. Usually 10m for chl",
    )
    parser.add_argument(
        "--var",
        "-v",
        type=str,
        required=True,
        choices=[
            "P_l",
            "kd490",
            "P1l",
            "P2l",
            "P3l",
            "P4l",
            "RRS412",
            "RRS443",
            "RRS490",
            "RRS510",
            "RRS555",
            "RRS670",
        ],
        help=""" model var name""",
    )
    parser.add_argument(
        "--datestart",
        "-t",
        type=date_from_str,
        required=True,
        help=""" Date start for time interval to consider for validation, format %Y%m%d""",
    )
    parser.add_argument(
        "--dateend",
        "-e",
        type=date_from_str,
        required=True,
        help=""" Date end for time interval to consider for validation,format %Y%m%d""",
    )
    parser.add_argument(
        "--zone",
        "-z",
        type=str,
        required=False,
        default="Med",
        help=""" Areas to generate the STATISTICS mean. std, bias and RMSD with respect satellite: Med or rivers""",
    )
    return parser.parse_args()


args = argument()

comm = get_mpi_communicator()
rank = comm.Get_rank()
nranks = comm.size


def weighted_mean_std(Conc, Weight):
    Weight_sum = Weight.sum()
    Mass = (Conc * Weight).sum()
    Weighted_Mean = Mass / Weight_sum
    Weighted_Std = np.sqrt(
        ((Conc - Weighted_Mean) ** 2 * Weight).sum() / Weight_sum
    )
    return Weighted_Mean, Weighted_Std


# def weighted_var(Conc, Weight):
#     ConcMean = weighted_mean(Conc, Weight)
#     Weight_sum = Weight.sum()
#     Mass = ((Conc - ConcMean) ** 2 * Weight).sum()
#     return Mass / Weight_sum

TheMask = Mask.from_file(args.maskfile)
Sup_mask = TheMask.cut_at_level(0)
MODEL_DIR = args.inputmodeldir
REF_DIR = args.satdir
area = args.zone
modvarname = args.var
satvarname = SAT_VARS[modvarname]

TI = TimeInterval.fromdatetimes(args.datestart, args.dateend)
dateformat = "%Y%m%d"

sat_TL = TimeList.fromfilenames(
    None, REF_DIR, "*.nc", prefix="", dateformat=dateformat
)
model_TL = TimeList.fromfilenames(TI, MODEL_DIR, "*.nc", filtervar=modvarname)

suffix = sat_TL.filelist[0].name[8:]

if area == "Med":
    from bitsea.basins import V2 as OGS
if area == "rivers":
    from bitsea.basins import RiverBoxes as OGS
if area == "coastal":
    from bitsea.basins import COASTAL12nm as OGS

BASINS = OGS.P

nFrames = model_TL.nTimes
nSub = len(OGS.P.basin_list)

jpk, jpj, jpi = TheMask.shape
dtype = [(sub.name, bool) for sub in OGS.P]
SUB = np.zeros((jpj, jpi), dtype=dtype)


for sub in BASINS:
    print(sub.name)
    sbmask = SubMask(sub, TheMask)
    SUB[sub.name] = sbmask[0, :, :]


mask200_2D = TheMask.mask_at_level(200.0)
mask0_2D = TheMask.mask_at_level(0.0)
# Do not consider Atlantic in the Mediterranean domain:
if area == "Med":
    SUB["med"] = mask0_2D.copy()
    ii = SUB["atl"]
    SUB["med"][ii] = False

COASTNESS_LIST = args.coastness
nCoast = len(COASTNESS_LIST)
COASTNESS = np.ones(
    (jpj, jpi), dtype=[(coast, bool) for coast in COASTNESS_LIST]
)
if "coast" in COASTNESS_LIST:
    COASTNESS["coast"] = coastmask = mask0_2D & (~mask200_2D)
if "open_sea" in COASTNESS_LIST:
    COASTNESS["open_sea"] = mask200_2D
if "everywhere" in COASTNESS_LIST:
    COASTNESS["everywhere"] = mask0_2D


# This is the surface layer choosen to match satellite chl data
surf_layer = Layer(0, args.layer)

modeltimes = model_TL.Timelist[rank::nranks]
filelist = model_TL.filelist[rank::nranks]
for itime, modeltime in enumerate(modeltimes):
    outfile = (
        args.outdir / f"valid.{modeltime.strftime(dateformat)}.{modvarname}.nc"
    )
    print(outfile, flush=True)
    modfile = filelist[itime]
    try:
        CoupledList = sat_TL.couple_with([modeltime])
        sattime = CoupledList[0][0]
    except IndexError:
        print("No sat file compliant with ", modfile)
        continue

    satfile = REF_DIR / (sattime.strftime(dateformat) + suffix)

    De = DataExtractor(TheMask, filename=modfile, varname=modvarname)
    Model = MapBuilder.get_layer_average(De, surf_layer)
    Sat = Satmodule.readfromfile(satfile, var=satvarname)
    cloudsLand = (np.isnan(Sat)) | (Sat > 1.0e19) | (Sat < 0)
    modelLand = np.isnan(Model)  # lands are nan
    nodata = cloudsLand | modelLand

    D = {
        "modelfile": str(modfile),
        "satfile": str(satfile),
        "COASTNESS_LIST": COASTNESS_LIST,
        "basinlist": [sub.name for sub in OGS.P],
    }
    D["VALID_POINTS"] = np.zeros((nSub, nCoast), np.float32)
    D["MODEL_MEAN"] = np.zeros((nSub, nCoast), np.float32)
    D["SAT___MEAN"] = np.zeros((nSub, nCoast), np.float32)
    D["MODEL__STD"] = np.zeros((nSub, nCoast), np.float32)
    D["SAT____STD"] = np.zeros((nSub, nCoast), np.float32)
    D["BGC_CLASS4_CHL_BIAS_SURF_BASIN"] = np.zeros((nSub, nCoast), np.float32)
    D["BGC_CLASS4_CHL_RMS_SURF_BASIN"] = np.zeros((nSub, nCoast), np.float32)
    D["BGC_CLASS4_CHL_CORR_SURF_BASIN"] = np.zeros((nSub, nCoast), np.float32)
    D["BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG"] = np.zeros(
        (nSub, nCoast), np.float32
    )
    D["BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG"] = np.zeros(
        (nSub, nCoast), np.float32
    )
    for icoast, coast in enumerate(COASTNESS_LIST):
        for isub, sub in enumerate(BASINS):
            selection = SUB[sub.name] & (~nodata) & COASTNESS[coast]
            M = matchup.matchup(Model[selection], Sat[selection])
            D["VALID_POINTS"][isub, icoast] = M.number()
            if M.number() == 0:
                continue
            D["BGC_CLASS4_CHL_RMS_SURF_BASIN"][isub, icoast] = M.RMSE()
            D["BGC_CLASS4_CHL_BIAS_SURF_BASIN"][isub, icoast] = M.bias()
            D["BGC_CLASS4_CHL_CORR_SURF_BASIN"][isub, icoast] = M.correlation()

            weight = TheMask.area[selection]
            D["MODEL_MEAN"][isub, icoast], D["MODEL__STD"][isub, icoast] = (
                weighted_mean_std(M.Model, weight)
            )
            D["SAT___MEAN"][isub, icoast], D["SAT____STD"][isub, icoast] = (
                weighted_mean_std(M.Ref, weight)
            )
            Mlog = matchup.matchup(
                np.log10(Model[selection]), np.log10(Sat[selection])
            )  # add matchup based on logarithm
            D["BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG"][isub, icoast] = Mlog.RMSE()
            D["BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG"][isub, icoast] = Mlog.bias()

    netcdf_validation_file.write(outfile, D)

import argparse

import netCDF4 as nc
import numpy as np
import numpy.ma as ma

from bitsea.basins import V2 as basV2
from bitsea.commons.layer import Layer
from bitsea.commons.mask import Mask
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.commons.utils import addsep
from bitsea.commons.utils import writetable
from bitsea.timeseries.plot import Hovmoeller_matrix
from bitsea.timeseries.plot import read_pickle_file


def argument():
    parser = argparse.ArgumentParser(
        description="""
    Generates in output directory two files ( model and ref) for each pfts (P1l, P2l, P3l, P4l)
    containing [nSub, nLayers] arrays of climatologies.

    if (P_type=="P1l"): P_name="DIATO"
    if (P_type=="P2l"): P_name="NANO"
    if (P_type=="P3l"): P_name="PICO"
    if (P_type=="P4l"): P_name="DINO"

    The files have the following nomenclature:
    infile_MEAN=INDIR + P_type + "_OPTION3_newSeasons3_means.csv"
    infile_STD=INDIR  + P_type + "_OPTION3_newSeasons3_std.csv"

    These arrays will be used in the next step to generate the following metrics:

    # Input: DIRECTORYof StatProfiles of P1l, P2l, P3l, P4l (files *.pkl)

    P1l-LAYER-Y-CLASS4-CLIM-BIAS/RMSD
    P2l-LAYER-Y-CLASS4-CLIM-BIAS/RMSD
    P3l-LAYER-Y-CLASS4-CLIM-BIAS/RMSD
    P4l-LAYER-Y-CLASS4-CLIM-BIAS/RMSD
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--inputdir",
        "-i",
        type=str,
        required=True,
        help="Directory of StatProfiles of P1l, P2l, P3l, P4l in *.pkl",
    )

    parser.add_argument(
        "--outdir", "-o", type=str, default=None, required=True, help=""
    )
    parser.add_argument(
        "--maskfile", "-m", type=str, default=None, required=True, help=""
    )
    parser.add_argument(
        "--starttime",
        "-s",
        type=str,
        required=True,
        help="start date in yyyymmdd format",
    )
    parser.add_argument(
        "--endtime",
        "-e",
        type=str,
        required=True,
        help="start date in yyyymmdd format",
    )
    return parser.parse_args()


args = argument()

LayerList = [
    Layer(0, 10),
    Layer(10, 30),
    Layer(30, 60),
    Layer(60, 100),
    Layer(100, 150),
    Layer(150, 300),
    Layer(0, 300),
]
rows_names = [layer.string() for layer in LayerList]

INPUTDIR = addsep(args.inputdir)
OUTDIR = addsep(args.outdir)
TI = TimeInterval(args.starttime, args.endtime, "%Y%m%d")

TheMask = Mask.from_file(args.maskfile)
jpk, jpj, jpi = TheMask.shape
z = -TheMask.zlevels

z_clim = np.array([-(el.bottom + el.top) / 2 for el in LayerList])

# DATASET of PFTS "CLIMATOLOGY" from HPLC
HPLC_CLIM_path = (
    "/g100_scratch/userexternal/lfeudale/HPLC/HPLC_pfts_CLIMATOLOGY/"
)
# HPLC_CLIM_path="/g100_scratch/userexternal/gocchipi/analisi_file_eva/"


def Layers_Mean(Pres, Values, LayerList):
    """
    Performs mean of profile along layers.

    Arguments:
    * Pres      * numpy array of pressures
    * Values    * numpy array of concetrations
    * LayerList * list of layer objects
    Returns :
    * MEAN_LAY * [nLayers] numpy array, mean of each layer
    """
    MEAN_LAY = np.zeros(len(LayerList), np.float32)
    STD_LAY = np.zeros(len(LayerList), np.float32)

    for ilayer, layer in enumerate(LayerList):
        ii = (Pres >= layer.top) & (Pres <= layer.bottom)
        if ii.sum() > 1:
            local_profile = Values[ii]
            MEAN_LAY[ilayer] = np.mean(local_profile)
            STD_LAY[ilayer] = np.std(local_profile)
    return MEAN_LAY, STD_LAY


#########

VARLIST = ["P1l", "P2l", "P3l", "P4l"]
SUBlist = basV2.Pred.basin_list
SUBlist.remove(SUBlist[-1])  # REMOVE ATLANTIC BUFFER
basin_names = [sub.name for sub in SUBlist]

nSub = len(SUBlist)
nLayers = len(LayerList)
METRICvar = {
    "P1l": "Diatoms",
    "P2l": "Nanophytoplankton",
    "P3l": "Picophytoplankton",
    "P4l": "Dinoflagellates",
}

VARnameNC = {
    "P1l": "media_seasonal_fdiatom_ALL",
    "P2l": "media_seasonal_fnano_ALL",
    "P3l": "media_seasonal_fpico_ALL",
    "P4l": "media_seasonal_fdino_ALL",
}


# media_seasonal_fdiatom, media_seasonal_fnano_ALL, media_seasonal_fpico_ALL, media_seasonal_fdino_ALL
# Remove Altantic Buffer from the list:

RESULTS_MEAN = np.zeros(
    (nLayers, 4), np.float32
)  # "meanWref","meanWmod","meanSref","meanSmod" --> 4 COLUMNS
RESULTS_STD = np.zeros(
    (nLayers, 4), np.float32
)  # "stdWref","stdWmod","stdSref","stdSmod"      --> 4 COLUMNS

rows_names = [layer.string() for layer in LayerList]


for ivar, var in enumerate(VARLIST):
    filename = INPUTDIR + var + ".pkl"
    TIMESERIES, TL = read_pickle_file(filename)

    print(METRICvar[var] + "-LAYER-Y-CLASS4-CLIM-BIAS,RMSD")

    BIAS = np.zeros((2, nLayers))
    RMSD = np.zeros((2, nLayers))
    CORR = np.zeros((2, nLayers))

    ########## DEFINE SEASONS: ##########
    for iSeas in [0, 1]:  # 0=Summer, 1=Winter
        if iSeas == 0:
            SeasonName = "summer"  # imonths in [5,6,7,8,9,10]
            ind_SUMMER = [
                ii
                for ii, x in enumerate(TL.Timelist)
                if x.month in [5, 6, 7, 8, 9, 10]
            ]
            ind_Seas = ind_SUMMER
        if iSeas == 1:
            SeasonName = "winter"  # imonths in [1,2,3,4,11,12]
            ind_WINTER = [
                ii
                for ii, x in enumerate(TL.Timelist)
                if x.month in [1, 2, 3, 4, 11, 12]
            ]
            ind_Seas = ind_WINTER

        TIMESERIES_season = TIMESERIES[
            ind_Seas, :, :, :, :
        ]  # dim(StatProf)==> [nFrames,nSub,nCoast,depth,nStat]
        TIMES = []

        for ii in ind_Seas:
            TIMES.append(TL.Timelist[ii])
        TLL = TimeList(TIMES)

        CLIM_MODEL = np.zeros((nSub, nLayers))
        CLIM_MODEL_STD = np.zeros((nSub, nLayers))

        Model_Med_Mean_Std = np.zeros((nLayers, 2), np.float32)

        for iSub, sub in enumerate(SUBlist):
            print(sub.name)
            Mean_profiles, _, _ = Hovmoeller_matrix(
                TIMESERIES_season, TLL, np.arange(jpk), iSub, icoast=1, istat=0
            )  # istat=0 --> MEAN
            mean_profile = Mean_profiles.mean(axis=1)
            mean_profile[mean_profile == 0] = np.nan

            Std_profiles, _, _ = Hovmoeller_matrix(
                TIMESERIES_season, TLL, np.arange(jpk), iSub, icoast=1, istat=1
            )  # istat=1 --> STD

            CLIM_MODEL[iSub, :], _ = Layers_Mean(
                TheMask.zlevels, mean_profile, LayerList
            )
            CLIM_MODEL_STD[iSub, :], _ = Layers_Mean(
                TheMask.zlevels, Std_profiles, LayerList
            )

        Model_Med_Mean_Std[:, 0] = np.nanmean(CLIM_MODEL, axis=0)
        Model_Med_Mean_Std[:, 1] = np.nanmean(CLIM_MODEL_STD, axis=0)
        writetable(
            OUTDIR + var + "_ModelMean_table_LayBas_" + SeasonName + ".txt",
            CLIM_MODEL,
            basin_names,
            rows_names,
        )
        writetable(
            OUTDIR + var + "_ModelStd_table_LayBas_" + SeasonName + ".txt",
            CLIM_MODEL_STD,
            basin_names,
            rows_names,
        )
        writetable(
            OUTDIR + var + "_ModelMean_table_Lay_" + SeasonName + ".txt",
            Model_Med_Mean_Std,
            rows_names,
            ["Mean", "Std"],
        )

        CLIM_REF_static = np.zeros(
            (2, nSub, nLayers, 2), np.float32
        )  # [SEASONS,NSUB,NLAYER,NMETRICS]

        # Read Climatology:
        infile = HPLC_CLIM_path + "OUTPUT" + VARnameNC[var] + ".nc"
        f = nc.Dataset(infile, "r")
        CLIM_REF_static[:, :, :, 0] = f.variables["mean_" + VARnameNC[var]][:]
        CLIM_REF_static[:, :, :, 1] = f.variables["std_" + VARnameNC[var]][:]
        f.close()

        Ref_Med_Mean_Std = np.zeros(
            (nLayers, 2), np.float32
        )  # LAYERS and SEASONS
        Ref_Med_Mean_Std[:, 0] = np.nanmean(
            CLIM_REF_static[iSeas, :, :, 0], axis=0
        )
        Ref_Med_Mean_Std[:, 1] = np.nanmean(
            CLIM_REF_static[iSeas, :, :, 1], axis=0
        )
        writetable(
            OUTDIR + var + "_CLIM_table_LayBas_" + str(iSeas) + ".txt",
            CLIM_REF_static[iSeas, :, :, 0],
            basin_names,
            rows_names,
        )  # ,["Mean","Std"])
        writetable(
            OUTDIR + var + "_CLIM_table_LayBas_sd" + str(iSeas) + ".txt",
            CLIM_REF_static[iSeas, :, :, 1],
            basin_names,
            rows_names,
        )
        writetable(
            OUTDIR + var + "_RefMean_table_Lay_" + SeasonName + ".txt",
            Ref_Med_Mean_Std,
            rows_names,
            ["Mean", "Std"],
        )

        BIAS[iSeas, :] = np.nanmean(
            CLIM_MODEL - CLIM_REF_static[iSeas, :, :, 0], axis=0
        )
        RMSD[iSeas, :] = np.sqrt(
            np.nanmean(
                (CLIM_MODEL - CLIM_REF_static[iSeas, :, :, 0]) ** 2, axis=0
            )
        )
        for ilayer, lay in enumerate(LayerList):
            CORR[iSeas, ilayer] = (
                ma.corrcoef(
                    ma.masked_invalid(CLIM_MODEL[:, ilayer]),
                    ma.masked_invalid(CLIM_REF_static[iSeas, :, ilayer, 0]),
                )
            )[0, 1]

        if iSeas == 1:  # WINTER
            print("WIN")
            RESULTS_MEAN[:, 0] = Ref_Med_Mean_Std[:, 0]
            RESULTS_MEAN[:, 1] = Model_Med_Mean_Std[:, 0]
            RESULTS_STD[:, 0] = Ref_Med_Mean_Std[:, 1]
            RESULTS_STD[:, 1] = Model_Med_Mean_Std[:, 1]

        if iSeas == 0:  # SUMMER
            print("SUM")
            RESULTS_MEAN[:, 2] = Ref_Med_Mean_Std[:, 0]
            RESULTS_MEAN[:, 3] = Model_Med_Mean_Std[:, 0]
            RESULTS_STD[:, 2] = Ref_Med_Mean_Std[:, 1]
            RESULTS_STD[:, 3] = Model_Med_Mean_Std[:, 1]

    writetable(
        OUTDIR + var + "_BIAS" + ".txt",
        np.transpose(BIAS),
        rows_names,
        ["summer", "winter"],
    )  # , basin_names,rows_names)
    writetable(
        OUTDIR + var + "_RMSD" + ".txt",
        np.transpose(RMSD),
        rows_names,
        ["summer", "winter"],
    )
    writetable(
        OUTDIR + var + "_CORR" + ".txt",
        np.transpose(CORR),
        rows_names,
        ["summer", "winter"],
    )

    writetable(
        OUTDIR + var + "_OPTION3_newSeasons3_means.csv",
        RESULTS_MEAN,
        rows_names,
        ["meanWref", "meanWmod", "meanSref", "meanSmod"],
    )
    writetable(
        OUTDIR + var + "_OPTION3_newSeasons3_std.csv",
        RESULTS_STD,
        rows_names,
        ["stdWref", "stdWmod", "stdSref", "stdSmod"],
    )

Cols_names = ["meanWref", "meanWmod", "meanSref", "meanSmod"]
Rows_names = [
    "(0,10)",
    "(10,30)",
    "(30,60)",
    "(100,150)",
    "(150,300)",
    "(0,300)",
]

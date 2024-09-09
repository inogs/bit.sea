import argparse

def argument():
    parser = argparse.ArgumentParser(description = 'Executes extraction of HC metrics vs SAT')

    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = 'Where the metrics HC are saved'
                                )

    parser.add_argument(   '--outputdir', '-o',
                                type = str,
                                required = True,
                                help = 'The directory where you want to dump compressed files'
                                )

    parser.add_argument(   '--version', '-v',
                                type = str,
                                required = True,
                                help = 'Name of the version of the actual chain (V9C, V10C)'
                                )


    return parser.parse_args()


import numpy as np
import netCDF4 as NC
from bitsea.commons.Timelist import TimeList
from bitsea.basins import V2 as OGS
from bitsea.commons.utils import writetable
from bitsea.commons.utils import addsep

args = argument()

INDIR = addsep(args.inputdir)
OUTDIR = addsep(args.outputdir)
version = args.version

SATtypes = ['DT','NRT']


# QuID 006-014 V1.3 (V6C)
EAN_RMS = {
    'win': 0.10,
    'sum': 0.06,
}
EAN_BIAS = {
    'win': 0.06,
    'sum': 0.01,
}

SeasonMonths = {
    'win': [11,12,1,2,3,4],
    'sum': [5,6,7,8,9,10],
}

nSub = len(OGS.P.basin_list)
for isub,sub in enumerate(OGS.P.basin_list):
    if 'med' in sub.name:
        index_med = isub

col_names=["BIAS","RMSE","EAN_bias","EAN_rmse"]

for tt in SATtypes:
    TL = TimeList.fromfilenames(None,INDIR,"Validation_hc_*" + tt + "*.nc", \
        prefix="Validation_hc_YYYYMMDD_on_weekly_Sat" + tt + ".", \
        dateformat='%Y%m%d')
    lenTL = TL.nTimes
    BIAS = np.zeros((lenTL,nSub))
    RMSE = np.zeros((lenTL,nSub))
    EAN_bias = np.zeros((lenTL,nSub))
    EAN_rmse = np.zeros((lenTL,nSub))
    Dates = []

    for ii,filein in enumerate(TL.filelist):
        dd = TL.Timelist[ii]
        datef = dd.strftime('%Y-%m-%d')
        #print(tt + ' ' + datef)
        Dates.append(datef)
        for ss in SeasonMonths.keys():
            if dd.month in SeasonMonths[ss]:
                seas = ss
        M = NC.Dataset(filein,"r")
        model = M.variables['MODEL_MEAN_LOG'][:,1]
        sat   = M.variables['SAT___MEAN_LOG'][:,1]
        BIAS[ii,:] = model-sat
        RMSE[ii,:] = M.variables['BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG'][:,1]
        EAN_bias[ii,:] = EAN_BIAS[seas]
        EAN_rmse[ii,:] = EAN_RMS[seas]

    startdate = TL.Timelist[0].strftime('%Y%m%d')
    enddate = TL.Timelist[-1].strftime('%Y%m%d')
    filetxt = OUTDIR + '/table_statistics_' + version + '_' + tt + '.txt'


    row_names=Dates
    METRICS=np.zeros((lenTL,4))*np.nan
    METRICS[:,0] = BIAS[:,index_med]
    METRICS[:,1] = RMSE[:,index_med]
    METRICS[:,2] = EAN_bias[:,index_med]
    METRICS[:,3] = EAN_rmse[:,index_med]

    writetable(filetxt,METRICS,row_names,col_names,fmt="%5.4f\t ") 

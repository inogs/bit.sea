# OUTPUTS
# images for QUID
# CMEMS-Med-QUID-006-008-V2-V1.0.docx
# Figure IV.3 and table IV.1

import argparse
from bitsea.utilities.argparse_types import date_from_str
from bitsea.utilities.argparse_types import existing_dir_path

def argument():
    parser = argparse.ArgumentParser(description = '''
    Plot timeseries fo BIAS and RMSE for QUID
    CMEMS-Med-QUID-006-008-V2-V1.0.docx
    Figure IV.3 and table IV.1
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = existing_dir_path,
                            required =True,
                            help = ''' Output image dir'''
                            )
    parser.add_argument(   '--inputdir', '-i',
                            type = existing_dir_path,
                            required = True,
                            help = 'Dir containing NetCDF validation files')
    parser.add_argument(   '--datestart', '-s',
                                type = date_from_str,
                                required =True,
                                help = ''' Date start for time interval to consider for validation, format %Y%m%d'''
                                )
    parser.add_argument(   '--dateend', '-e',
                                type = date_from_str,
                                required =True,
                                help = ''' Date end for time interval to consider for validation,format %Y%m%d'''
                                )
    parser.add_argument(   '--var', '-v',
                                type = str,
                                required = True,
                                choices = ['P_l','kd490','P1l','P2l','P3l','P4l','RRS412','RRS443','RRS490','RRS510','RRS555','RRS670' ],
                                help = ''' model var name'''
                                )
    parser.add_argument(   '--coastness', '-c',
                                type = str,
                                choices=["coast","open_sea","everywhere"]
                                )
    parser.add_argument(   '--zone', '-z',
                                type = str,
                                required =False,
                                default = "Med",
                                help = ''' Areas to generate the STATISTICS mean. std, bias and RMSD with respect satellite: Med or rivers'''
                                )



    return parser.parse_args()

args = argument()


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import numpy as np
from bitsea.instruments.var_conversions import SAT_VARS
from bitsea.commons.Timelist import TimeList, TimeInterval
from bitsea.validation.deliverables import netcdf_validation_file



if (args.zone == "Med"):
    from bitsea.basins import V2 as OGS
if (args.zone == "rivers"):
    print ("rivers")
    from bitsea.basins import RiverBoxes as OGS

model_label=' MODEL'

if (args.var =="kd490"):
    units="[m$^{-1}$]"
    vmin=-0.1
    vmax=0.1
elif (args.var.startswith('RRS')):
    units="[st$^{-1}$]"
    vmin=-0.01 
    vmax=0.01 
else:
    units="[mg/m$^3$]"
    vmin=-0.3
    vmax=0.3
    if (args.zone == "rivers"):
        vmin=-1.0
        vmax=1.0
    
var_label = SAT_VARS[args.var] + " " + units

TI = TimeInterval.fromdatetimes(args.datestart, args.dateend)
dr=netcdf_validation_file.dir_reader(TI,args.inputdir,args.var,args.coastness)



nSUB = len(OGS.P.basin_list)

for isub,sub in enumerate(OGS.P):
    if (sub.name == 'atl') : continue
    outfilename=args.outdir / f"{args.var}-RMS-BIAS_{sub.name}.png"
    print (outfilename)
    if ((args.var =="P_l") & (sub.name=="Po")):
        vmin=-4.0
        vmax=4.0

    fig, ax = pl.subplots()
    fig.set_size_inches(12,4)
    ax.plot(dr.TIMES,dr.BGC_CLASS4_CHL_RMS_SURF_BASIN[:,isub],'-k',label='RMS')
    ax.plot(dr.TIMES,dr.BGC_CLASS4_CHL_BIAS_SURF_BASIN[:,isub],'-b',label='Bias')
    ax.set_ylabel(sub.name.upper() + ' - ' + var_label).set_fontsize(14)
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    pl.setp(ltext,fontsize=12)
    pl.rc('xtick', labelsize=12)
    pl.rc('ytick', labelsize=12)
    pl.ylim(vmin,vmax)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
    ax.grid(True)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30,fontsize=9)
    ylabels = ax.get_yticklabels()
    pl.setp(ylabels, fontsize=10)
    pl.tight_layout()
    pl.tight_layout()
    pl.savefig(outfilename)
    pl.close(fig)


from bitsea.commons.season import season
S=season()
S.setseasons(["0101", "0501", "0601", "1001"], ["winter","spring","summer","fall"])
from bitsea.commons import timerequestors
TL=TimeList(dr.TIMES)
from bitsea.commons.utils import writetable

iSeas=0 # JAN-APR
CLIM_REQ=timerequestors.Clim_season(iSeas,S)
ii,w=TL.select(CLIM_REQ)
RMS__win = np.nanmean(dr.BGC_CLASS4_CHL_RMS_SURF_BASIN[     ii,:],axis=0)
BIAS_win = np.nanmean(dr.BGC_CLASS4_CHL_BIAS_SURF_BASIN[    ii,:],axis=0)
RMSL_win = np.nanmean(dr.BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[ ii,:],axis=0)
BIASLwin = np.nanmean(dr.BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[ii,:],axis=0)

MEAN_MOD_win = np.nanmean(dr.MODEL_MEAN[ii,:],axis=0)
MEAN_REF_win = np.nanmean(dr.SAT___MEAN[ii,:],axis=0)

STD_MOD_win = np.nanmean(dr.MODEL_STD[ii,:],axis=0)
STD_REF_win = np.nanmean(dr.SAT___STD[ii,:],axis=0)
CORR_win    = np.nanmean(dr.BGC_CLASS4_CHL_CORR_SURF_BASIN[ii,:],axis=0)

iSeas=2 # JUN-SEP
CLIM_REQ=timerequestors.Clim_season(iSeas,S)
ii,w=TL.select(CLIM_REQ)
RMS__sum = np.nanmean(dr.BGC_CLASS4_CHL_RMS_SURF_BASIN[     ii,:],axis=0)
BIAS_sum = np.nanmean(dr.BGC_CLASS4_CHL_BIAS_SURF_BASIN[    ii,:],axis=0)
RMSL_sum = np.nanmean(dr.BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[ ii,:],axis=0)
BIASLsum = np.nanmean(dr.BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[ii,:],axis=0)

MEAN_MOD_sum = np.nanmean(dr.MODEL_MEAN[ii,:],axis=0)
MEAN_REF_sum = np.nanmean(dr.SAT___MEAN[ii,:],axis=0)

STD_MOD_sum = np.nanmean(dr.MODEL_STD[ii,:],axis=0)
STD_REF_sum = np.nanmean(dr.SAT___STD[ii,:],axis=0)
CORR_sum    = np.nanmean(dr.BGC_CLASS4_CHL_CORR_SURF_BASIN[ii,:],axis=0)

mat = np.zeros((nSUB,18),np.float32)

mat[:,0] = RMS__win
mat[:,1] = RMS__sum
mat[:,2] = BIAS_win
mat[:,3] = BIAS_sum
mat[:,4] = RMSL_win
mat[:,5] = RMSL_sum
mat[:,6] = BIASLwin
mat[:,7] = BIASLsum
mat[:,8] = STD_MOD_win
mat[:,9] = STD_REF_win
mat[:,10] = STD_MOD_sum
mat[:,11] = STD_REF_sum
mat[:,12] = CORR_win
mat[:,13] = CORR_sum
#----
mat[:,14] = MEAN_MOD_win
mat[:,15] = MEAN_REF_win
mat[:,16] = MEAN_MOD_sum
mat[:,17] = MEAN_REF_sum

outfiletable = args.outdir / ("table4.1_" + args.var + ".txt")
print (outfiletable)
rows_names=[sub.name for sub in OGS.P.basin_list]
column_names=['RMSwin','RMSsum','BIASwin','BIASsum', 'RMSLwin','RMSLsum','BIASLwin','BIASLsum','STD_MODwin','STD_SATwin','STD_MODsum','STD_SATsum','CORRwin','CORRsum','MEAN_MODwin','MEAN_SATwin','MEAN_MODsum','MEAN_SATsum']
writetable(outfiletable, mat, rows_names, column_names, fmt='%5.3f\t')


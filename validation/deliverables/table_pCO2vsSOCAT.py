import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Produces tables [sub, month] of means of surface values in 2017
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = False,
                                help = '''DIRECTORY of MONTHLY MODEL DATA''')

#    parser.add_argument(   '--outdir', '-o',
#                                type = str,
#                                default = None,
#                                required = True,
#                                help = "")

    return parser.parse_args()

args = argument()

from commons.utils import addsep
INPUTDIR=addsep(args.inputdir)
filemodel=INPUTDIR + "monthly_pCO2.txt"
import numpy as np
from basins import V2
from commons.utils import writetable

socat=np.loadtxt("monthly_clim_socat.txt",skiprows=1,usecols=range(1,13))
num_socat=np.loadtxt("monthly_num_socat.txt",skiprows=1,usecols=range(1,13))
model=np.loadtxt(filemodel,skiprows=1,usecols=range(1,13))

nSUB = len(V2.P.basin_list)
M = np.zeros((nSUB, 4),np.float32)*np.nan
TOT_RMSD = np.zeros((1, 2),np.float32)*np.nan

for isub, sub in enumerate(V2.P.basin_list):
    BIAS=np.nanmean(model[isub,:]-socat[isub,:])
    RMSD=np.sqrt(np.nanmean((model[isub,:]-socat[isub,:])**2))
    bad = ~np.logical_or(np.isnan(model[isub,:]), np.isnan(socat[isub,:]))
    mod1=np.compress(bad, model[isub,:])  # array([  5.,   1.,   6.,  10.,   1.,   1.])
    soc1=np.compress(bad, socat[isub,:])
    C=np.corrcoef(mod1,soc1)
    n_socat=np.nansum(num_socat[isub,:])

    M[isub,0]=BIAS
    M[isub,1]=RMSD
    M[isub,2]=C[0,1]
    M[isub,3]=n_socat
    print sub.name, " ","%8.2f"%  BIAS, " ","%8.2f"%  RMSD, " ", "%8.2f"%  C[0,1], " ", "%d"%  n_socat
    # C=np.corrcoef(model[iSub,:],socat[iSub,:])

TOT_RMSD[0,0]=np.sqrt(np.nanmean((model[:,:-1]-socat[:,:-1])**2))
TOT_RMSD[0,1]=np.nansum(num_socat[:,:-1])

print "TOTAL RMSD: ", "%8.2f"%  TOT_RMSD[0,0]
rows_names_list=[sub.name for sub in V2.P]
column_names_list=["BIAS","RMSD","CORR","nSOCAT"]
outfile="pCO2-SURF-M-CLASS4-CLIM-RMSD-BASIN.txt"
writetable(outfile, M, rows_names_list, column_names_list,fmt=("%8.2f\t   %8.2f\t   %8.2f\t   %d"))
writetable("TOT_RMSD_pCO2vsSOCAT.txt",TOT_RMSD,["MED DOMAIN"],["TOT RMSD","NUM"],fmt=("%8.2f\t %d"))

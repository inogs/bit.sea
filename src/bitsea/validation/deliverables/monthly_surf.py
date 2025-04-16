import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Produces tables [sub, month] of means of surface values for a given year
    in files called monthl.var.txt
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = False,
                                help = '''A STAT_PROFILES dir with pkl files''')
    parser.add_argument(   '--year', '-y',
                                type = str,
                                required = True
                                )    

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

from bitsea.timeseries.plot import read_pickle_file
from bitsea.commons import timerequestors
from bitsea.basins import V2 as OGS
import numpy as np
from bitsea.commons.utils import writetable,addsep

INPUTDIR=addsep(args.inputdir)
OUTPUTDIR=addsep(args.outdir)
VARLIST=[ 'N1p', 'N5s', 'N3n', 'pH', 'ALK', 'DIC','pCO2' ]
#VARLIST=[ 'vosaline','votemper' ]

for var in VARLIST:
    filename=INPUTDIR + var + ".pkl"

# TL is determined by the filelists of STAT PROFILES in that directory (therefore on "file".pkl)
    data, TL = read_pickle_file(filename)
    nSUB = len(OGS.P.basin_list)

    MONTHLY = np.zeros((nSUB, 12),np.float32)*np.nan

    for imonth in range(12):
#        req=timerequestors.Monthly_req(int(args.year),imonth+1)
# Consider the month climatology:
        req=timerequestors.Clim_month(imonth+1)
        ii, w = TL.select(req)
        print (ii)
        print (w)
        for kk in ii: print (TL.Timelist[kk])
        for isub, sub in enumerate(OGS.P):
            V=data[ii,isub,1,0,0]
# They are weekly values, therefore we include the weight in the mean
#            MONTHLY[isub,imonth] = V.mean()
            MONTHLY[isub,imonth] = np.sum(w*V)/np.sum(w)

    rows_names_list=[sub.name for sub in OGS.P]
    column_names_list=[str(i) for i in range(1,13)]
    outfile="%smonthly_%s.txt" %(OUTPUTDIR,var)
    writetable(outfile, MONTHLY, rows_names_list, column_names_list)



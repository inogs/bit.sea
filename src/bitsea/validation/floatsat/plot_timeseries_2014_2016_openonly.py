# OUTPUTS
# images for QUID
# CMEMS-Med-QUID-006-008-V2-V1.0.docx
# Figure IV.2

import pickle
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import sys
import numpy as np
import argparse
#
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.mask import Mask
from commons.submask import SubMask
from commons.layer import Layer
from timeseries.plot import Hovmoeller_matrix
import numpy as np
from basins import V2
import matplotlib.pyplot as pl
from commons.utils import getcolor, addsep
#
def argument():
    parser = argparse.ArgumentParser(description = '''
    plot somethings
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    #parser.add_argument(   '--inputmodeldir', '-i',
    #                            type = str,
    #                            required =True,
    #                            help = ''' Input model dir, where P_l files are, usually ../wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/'''
    #                            )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')

    parser.add_argument(   '--outdir', '-O',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Output image dir'''
                            )
    parser.add_argument(   '--open_sea_file', '-o',
                            type = str,
                            required = True,
                            help = 'Input pickle file for open sea')

    


    return parser.parse_args()

args = argument()
#INPUTDIR = addsep(args.inputmodeldir)
fid = open(args.open_sea_file)
LIST = pickle.load(fid)
fid.close()
TIMES,_,_,MODEL_MEAN,SAT___MEAN,_,_ = LIST
    

TheMask = Mask(args.maskfile)
jpk,jpj,jpi=TheMask.shape

#TI = TimeInterval("201403","201501","%Y%m")
#TL =TimeList.fromfilenames(TI, INPUTDIR, "ave*nc")
var = 'P_l'

from basins import V2 as OGS
for isub,sub in enumerate(OGS.P):
    print sub.name
    fig, ax = pl.subplots()
#
    overall_isub = (OGS.P.basin_list).index(sub)
#    Mean_profiles,_,_ = Hovmoeller_matrix(TL.Timelist,TL.filelist, var, overall_isub, coast=0, stat=0, depths=np.arange(jpk)) #72 nFiles
#    fig,ax = pl.subplots()
#    ax.plot(TL.Timelist,Mean_profiles[0,:],'-k')
#
    ax.plot(TIMES,SAT___MEAN[:,isub],'og',label=' SAT')
    ax.plot(TIMES,MODEL_MEAN[:,isub],'-k',label=' MODEL open sea')
    ax.set_ylabel(sub.name.upper() + ' - CHL [mg/m$^3$]').set_fontsize(14)
    ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    pl.setp(ltext,fontsize=12)
    pl.rc('xtick', labelsize=12)
    pl.rc('ytick', labelsize=12)
    pl.ylim(0.0, 0.6)
    #ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
    ax.grid(True)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30)
    #ax.tick_params(direction='left', pad=2)
    #fig.show()
    outfilename=args.outdir+"/"+'chl' + sub.name + ".png"
    pl.savefig(outfilename)

import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    
    Creates Hovmoeller images, a png file for each subbasin for a selected variable.

    Caveats: here subbasin list here is hardcoded and must be consistent with the one used
    to generate STAT_PROFILES
    
    
    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = '/some/path/with/STAT_PROFILES')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = 'Output directory'
                                )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file'
                                )
    parser.add_argument(   '--varname', '-v',
                                type = str,
                                required = True,
                                help = 'name of the ogstm variable'
                                )
    return parser.parse_args()

args = argument()

from timeseries.plot import Hovmoeller_matrix, Hovmoeller_diagram
from timeseries.plot import read_pickle_file, read_basic_info
from basins import V2 as OGS
from commons.mask import Mask
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
from commons.utils import addsep


INPUTDIR=addsep(args.inputdir)
OUTDIR  =addsep(args.outdir)
TheMask = Mask(args.maskfile)

var=args.varname
filename=INPUTDIR + var + ".pkl"
iLev = TheMask.getDepthIndex(400) #m

TIMESERIES,TL=read_pickle_file(filename)
SUBLIST, COASTLIST, STAT_LIST = read_basic_info(TL.filelist[0])

icoast = COASTLIST.index('open_sea')
istat =  STAT_LIST.index('Mean')

years = mdates.YearLocator(1) #ticks every 10 years
months = mdates.MonthLocator()
yearsFmt = mdates.DateFormatter('%Y')



for isub, sub in enumerate(OGS.P):
    outfile = OUTDIR + 'Hovm_' + sub.name + '_' + var + '.png'
    print outfile
    M,xs,ys =  Hovmoeller_matrix(TIMESERIES, TL, TheMask.zlevels[:iLev+1], isub, icoast, istat)
    fig, ax0, im = Hovmoeller_diagram(M, xs, ys,shading="flat")# vmin=0, vmax=0.5)
    fig.colorbar(im)
    #ax0.xaxis.set_major_locator(years)
    
    #y1=datetime.datetime(1970,1,1).toordinal()
    #y2=datetime.datetime(2100,1,1).toordinal()
    #ax0.set_xlim([y1,y2])
    ax0.set_ylim([0,400])
    ax0.invert_yaxis()
    ax0.tick_params(axis='both', which='major', labelsize=15) 
    fig.suptitle(sub.name + ': ' + var)
#    fig.autofmt_xdate()
    fig.set_size_inches(30,8)
    fig.savefig(outfile)
    pl.close(fig)

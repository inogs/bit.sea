
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.mask import Mask
from plot import plot_Hovmoeller_diagram
from basins import V2 as OGS
import scipy.io.netcdf as NC


def read_basic_info(filename):
    ncIN = NC.netcdf_file(filename,'r')

    SUBLIST = ncIN.sub___list.split(", ")
    COASTLIST=ncIN.coast_list.split(",")
    STAT_LIST=ncIN.stat__list.split(",")    
    ncIN.close()
    return SUBLIST, COASTLIST, STAT_LIST

TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc')

INPUTDIR="/pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/POSTPROC/output/AVE_FREQ_2/PUB_SS/STAT_PROFILES/"

TI = TimeInterval('1990101','20150101',"%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*.nc")

SUBLIST, COASTLIST, STAT_LIST = read_basic_info(TL.filelist[0])


varname='P_i'
subbasin = 12
stat = 0
depths = 30
coast = 0
#fig, ax, im = plot_Hovmoeller_diagram(TL.filelist, varname, subbasin, coast, stat, TheMask.zlevels[:30])

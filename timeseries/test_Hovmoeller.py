
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.mask import Mask
from timeseries.plot import *
from basins import V2 as OGS
import scipy.io.netcdf as NC


def read_basic_info(filename):
    ncIN = NC.netcdf_file(filename,'r')

    SUBLIST = ncIN.sub___list.split(", ")
    COASTLIST=ncIN.coast_list.split(",")
    STAT_LIST=ncIN.stat__list.split(",")    
    ncIN.close()
    return SUBLIST, COASTLIST, STAT_LIST

TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc')

INPUTDIR="/pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/POSTPROC/output/AVE_FREQ_2/PUB_SS/STAT_PROFILES/"
#INPUTDIR="/pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/POSTPROC/output/AVE_FREQ_1/PUB_SS/STAT_PROFILES/"
TI = TimeInterval('1990101','20150101',"%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*.nc")

SUBLIST, COASTLIST, STAT_LIST = read_basic_info(TL.filelist[0])


varname='P_i'
subbasin = SUBLIST.index('tyr2')
stat = STAT_LIST.index('Mean')
coast = COASTLIST.index('open_sea')


#fig, ax, im = plot_Hovmoeller_diagram(TL.filelist, varname, subbasin, coast, stat, TheMask.zlevels[:30])

M,xs, ys = Hovmoeller_matrix(TL.filelist, varname, subbasin, coast, stat, TheMask.zlevels[:30])
fig, ax, im = Hovmoeller_diagram(M, xs, ys)
fig.colorbar(im)
fig.suptitle('Chlorophyll')
ax.set_ylabel('depth (m)')
#fig.autofmt_xdate()

fig.show()

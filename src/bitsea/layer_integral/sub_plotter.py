import numpy as np
from bitsea.basins import V2 as OGS
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
import matplotlib.pyplot as pl
from bitsea.layer_integral.mapplot import generic_mapplot_medeaf
import matplotlib.font_manager as font_manager
from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap

# generic_mapplot_medeaf needs
# ln -s /g100_work/OGS_prod2528/OPA/V12C-prod/etc/static-data/POSTPROC/fonts/TitilliumWeb-Regular.ttf
# ln -s /g100_work/OGS_prod2528/OPA/V12C-prod/etc/static-data/POSTPROC/fonts/TitilliumWeb-Bold.ttf


TheMask_all = Mask.from_file('/g100_work/OGS_prod2528/OPA/V12C-prod/etc/static-data/MED24_125/meshmask.nc')
TheMask = Mask(TheMask_all.grid, np.zeros((1,)), TheMask_all.mask[0,:,:],allow_broadcast=True)
jpk,jpj,jpi = TheMask.shape


SUB_matrix=np.zeros((jpj,jpi),np.float32)

SUBLIST = OGS.Pred.basin_list
nSub = len(SUBLIST)
xC = np.zeros(nSub,np.float32)
yC = np.zeros(nSub,np.float32)
colortext=['w','w','w','w','w','k','k','k','k','k','k','k','k','k','w','w','w']
cmap=pl.get_cmap('jet',nSub)

for isub, sub in enumerate(SUBLIST):
    m = SubMask(sub, TheMask)
    submask = m.mask[0,:,:]
    SUB_matrix[submask] =isub
    xC[isub]=TheMask.xlevels[submask].mean()
    yC[isub]=TheMask.ylevels[submask].mean()

SUB_matrix[~TheMask.mask[0,:,:]] = np.nan


xlim = [-6.5, 36.5]
ylim = [30, 46]
xc = (xlim[0] + xlim[1]) / 2
yc = (ylim[0] + ylim[1]) / 2

mapobj = Basemap(
    projection="merc",
    lat_0=xc,
    lon_0=yc,
    llcrnrlon=xlim[0],
    llcrnrlat=ylim[0],
    urcrnrlon=xlim[1],
    urcrnrlat=ylim[1],
    area_thresh=None,
    resolution="i",
)

map_dict = {'data':SUB_matrix, 'clim':[0,nSub-1]}

logo = pl.imread("/g100_work/OGS_prod2528/OPA/V12C-prod/etc/static-data/POSTPROC/LogoEchoOGS4.png")
pl.close('all')
fig,ax = generic_mapplot_medeaf(map_dict, mapobj, TheMask, fig=None, ax=None, ncolors=nSub, logo=logo)


font=FontProperties()
font_prop   = font_manager.FontProperties(fname='TitilliumWeb-Regular.ttf', size=12)
for isub,sub in enumerate(SUBLIST):
    x,y=mapobj(xC[isub],yC[isub])
    ax.annotate(sub.name, xy=(x,y), xytext=(0,0),textcoords="offset points", fontproperties=font_prop,color=colortext[isub], ha='center', va='center', alpha=1.0)


fig.savefig('Subbasin_V3_200_color.png',dpi=96)


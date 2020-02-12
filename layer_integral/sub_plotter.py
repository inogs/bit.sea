import numpy as np
from basins import V2 as OGS
from commons.mask import Mask
from commons.submask import SubMask
import matplotlib.pyplot as pl
from layer_integral.mapplot import generic_mapplot_medeaf
import matplotlib.font_manager as font_manager
from matplotlib.font_manager import FontProperties


TheMask_all=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc')
TheMask = TheMask_all.cut_at_level(0)
mask = TheMask.mask_at_level(0)
jpk,jpj,jpi = TheMask.shape


SUB_matrix=np.zeros((jpj,jpi),np.float32)

SUBLIST = OGS.Pred.basin_list
nSub = len(SUBLIST)
xC = np.zeros(nSub,np.float32)
yC = np.zeros(nSub,np.float32)
colortext=['w','w','w','w','w','k','k','k','k','k','k','k','k','k','w','w']
cmap=pl.get_cmap('jet',nSub)

for isub, sub in enumerate(SUBLIST):
    m = SubMask(sub,maskobject=TheMask)
    submask = m.mask[0,:,:]
    SUB_matrix[submask] =isub
    xC[isub]=TheMask.xlevels[submask].mean()
    yC[isub]=TheMask.ylevels[submask].mean()

SUB_matrix[~mask] = np.NaN

x,y = TheMask_all.coastline(200,50) # executed before axis creation, in order to avoid the plot of all contours
map_dict = {'data':SUB_matrix, 'clim':[0,nSub-1]}
bkg=pl.imread("/pico/home/usera07ogs/a07ogs00/OPA/V2C-dev/etc/static-data/POSTPROC/background_medeaf.png")
fig,ax = generic_mapplot_medeaf(map_dict, fig=None, ax=None, mask=TheMask, ncolors=nSub, background_img=bkg)
ax.plot(x,y,color='k',linewidth=0.7)

#pl.close('all')
font=FontProperties()
font_prop   = font_manager.FontProperties(fname='TitilliumWeb-Regular.ttf', size=12)
for isub,sub in enumerate(SUBLIST):
    t = ax.text(xC[isub],yC[isub],sub.name, color=colortext[isub], verticalalignment='center', horizontalalignment='center')
    t.set_font_properties(font_prop)

fig.savefig('Subbasin_V3_200_color.png',dpi=96)


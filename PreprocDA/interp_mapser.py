import numpy as np
from commons.interpolators import surf_interp_2d
from commons.mask import Mask

infile = 'mapser.19771215.16.npy'
outfile = 'mapser.19771215.4.npy'

# Loading masks
Mask16 = Mask('/pico/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc')
Mask4 = Mask('/pico/scratch/userexternal/ateruzzi/MASKS4/meshmask.nc')
#Mask24 = Mask('/pico/scratch/userexternal/ateruzzi/MASKS24/meshmask.nc', \
#              dzvarname='e3t_0')
maskcoast4 = (~Mask4.mask_at_level(200)) & Mask4.mask[0,:,:]


# reading mapser 1/16
mapser16 = np.load(infile)
mapser16[mapser16==0] = np.nan

# interpolation
mapser4 = surf_interp_2d(Mask16,Mask4,mapser16)

mapser4[~maskcoast4] = np.nan
mapser4[np.isnan(mapser4)] = 0

np.save(outfile,mapser4)

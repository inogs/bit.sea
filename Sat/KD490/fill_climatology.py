import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Fill climatology with nearest
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(   '--inputclima', '-i',
                                type = str,
                                required =True,
                                help = ''' Climatology input .nc'''
                                )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = 'Output dir ')
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')



    return parser.parse_args()

args = argument()


from commons.mask import Mask
from commons.utils import addsep
import numpy as np
import netCDF4
import scipy.interpolate

#TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc')
TheMask = Mask(args.maskfile)
#CLIM_FILE="Climatology_KD490.nc"
CLIM_FILE = args.inputclima
OUTDIR = addsep(args.outdir)

tmask = TheMask.mask_at_level(0)


ncIN = netCDF4.Dataset(CLIM_FILE,'r')
CLIMNN = np.array( ncIN.variables['Mean'])
ncIN.close()

jpk,jpj,jpi = TheMask.shape
#start_i=53

CLIMNN[CLIMNN==-999] = np.nan
#tmask=tmask[:, start_i:]
CLIMNN_FILLED=np.ones_like(CLIMNN)*-999




for julian in range(365):
    print julian
    #climNN=CLIMNN[julian,:,start_i:]
    climNN=CLIMNN[julian,:,:]
    goods  = (~np.isnan(climNN) & tmask)
    Jgoods, Igoods = np.nonzero(goods)
    nP = len(Jgoods)
    points = np.zeros((nP,2),dtype=np.float32)
    points[:,0] = Jgoods
    points[:,1] = Igoods
    values = climNN[goods]
    
    tofill = (np.isnan(climNN) & tmask)
    J,I = np.nonzero(tofill)
    nP = len(J)
    xi = np.zeros((nP,2),dtype=np.float32)
    xi[:,0] = J
    xi[:,1] = I
    
    V= scipy.interpolate.griddata(points, values, xi, "nearest")
    climNN[tofill] = V
    #CLIMNN_FILLED[julian,:,start_i:] = climNN
    CLIMNN_FILLED[julian,:,:] = climNN
    CLIMNN_FILLED[julian,tmask==0] = np.nan

CLIMNN_FILLED[np.isnan(CLIMNN_FILLED)] = -999.0

ncOUT = netCDF4.Dataset(OUTDIR + '/KD490_Climatology_24_filled.nc','w')
ncOUT.createDimension('lon',jpi)
ncOUT.createDimension('lat',jpj)
ncOUT.createDimension('time',365)
ncvar=ncOUT.createVariable('Mean','f',('time','lat','lon'))
setattr(ncvar,'missing_value',-999.0)
ncvar[:]=CLIMNN_FILLED
ncOUT.close()

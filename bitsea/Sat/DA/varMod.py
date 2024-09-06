from commons import netcdf4
from Sat import SatManager
from postproc import masks
from commons.mask import Mask
import numpy as np
import os,sys

VARERRDIR="/pico/scratch/userexternal/gbolzon0/EOF-python/VarErr/"
VARSATDIR="/pico/scratch/userexternal/gbolzon0/EOF-python/VARSAT/"

OUTPUTDIR="/pico/scratch/userexternal/gbolzon0/EOF-python/VarMod/"

maskfile="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/eas_v12/ogstm/meshmask.nc"

TheMask=Mask(maskfile)

minMod=0.5
minModC=0.75

mask_200= TheMask.mask_at_level(200)
mask_0  = TheMask.mask_at_level(0)

if (minModC == 0.5  ) : varversion="01"
if (minModC == 0.75 ) : varversion="02"

fillvalue = -999.0

for month in range(1,13):
    fileVE = "%svarErr.%02d.nc"  %( VARERRDIR, month)
    fileVS = "%svar2Dsat.CCI.F10.2.%02d.nc"  %( VARSATDIR, month)
    outfile = "%svar2Dcoast.%02dv%s.nc" %( OUTPUTDIR, month, varversion)
    VE = SatManager.readfromfile(fileVE,'variance')
    VS = SatManager.readfromfile(fileVS,'variance')
    VE[VE==fillvalue] = np.nan
    VS[VS==fillvalue] = np.nan
    varMod = VE - VS

    minMask_Open = ( varMod < VE*minMod ) & mask_200
    varMod[minMask_Open] = VE[minMask_Open]*minMod
    
    minMask_Coast = ( varMod < VE*minModC) & (~mask_200)
    varMod[minMask_Coast] = VE[minMask_Coast] * minModC


    JN,IN = np.where( np.isnan(varMod) & mask_0 )
    nP = len(JN)
    for ip in range(nP):
        jj = JN[ip]
        ji = IN[ip]
        varMod[jj, ji ] = np.nanmean(varMod[jj-1:jj+2, ji-1:ji+2] )
            
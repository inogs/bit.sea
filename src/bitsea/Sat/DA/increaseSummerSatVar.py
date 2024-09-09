from commons.mask import Mask
from Sat import SatManager
from postproc import masks
maskSat = getattr(masks,"Mesh24")

TheMask=Mask('/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/eas_v12/ogstm/meshmask.nc')
mask_200= TheMask.mask_at_level(200)

INPUTDIR="/pico/scratch/userexternal/gbolzon0/EOF-python/VARSAT/"
OUTPUTDIR="/pico/scratch/userexternal/gbolzon0/EOF-python/SUMMER_INCREASED50_100/"

incr = [4,16,16,16,16,4];  # April- September
incrCoast = [1,4,4,4,4,1];

for imonth, month in enumerate(range(4,10)) : # April- September
    inputfile ="%svar2Dsat.CCI.F10.2.%02d.nc"  %(  INPUTDIR, month)
    outfile   ="%svarSatinc%02d.nc"  %( OUTPUTDIR, month)
    varSat = SatManager.readfromfile(inputfile, "variance")
    ii = varSat == -999
    varSat[ mask_200] = varSat[ mask_200]*incr[imonth]
    varSat[~mask_200] = varSat[~mask_200]*incrCoast[imonth]
    varSat[ii] = -999
    SatManager.dumpGenericNativefile(outfile, varSat, "variance", maskSat)
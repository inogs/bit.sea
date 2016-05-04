from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import commons.IOnames as IOnames
import numpy as np
import SatManager as Sat
import matchup.matchup as matchup
import scipy.io.netcdf as NC


MODEL_DIR="/gpfs/work/OGS_prod/CalVal/Q_REP_MODEL_FORECAST/"
REF_DIR  = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/WEEKLY/"

Timestart="20150701"
Time__end="20151001"
TI    = TimeInterval(Timestart,Time__end,"%Y%m%d")
model_TL = TimeList.fromfilenames(TI, MODEL_DIR,"*.nc",prefix='',dateformat='%Y%m%d')

IonamesFile = '../postproc/IOnames_sat.xml'
IOname = IOnames.IOnames(IonamesFile)

ngib=52

for itime, time in enumerate(model_TL.Timelist[:1]):
    satfile = REF_DIR + time.strftime(IOname.Input.dateformat) + IOname.Output.suffix + ".nc"
    modfile = model_TL.filelist[itime]
     
    ncIN = NC.netcdf_file(modfile,'r')
    For = ncIN.variables['chl'].data[0,0,:,:].copy().astype(np.float64)
    ncIN.close()
    
    Sat16 = Sat.convertinV4format(Sat.readfromfile(satfile)).astype(np.float64)
    Sat16 = Sat16[:,ngib:]
    
    cloudsLand = np.isnan(Sat16)
    modelLand  = For > 1.0e+19
    nodata     = cloudsLand | modelLand
    M = matchup.matchup(For[~nodata], Sat16[~nodata])
    print M.bias()
    
    
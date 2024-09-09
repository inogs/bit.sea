from commons.timeseries import TimeSeries
from commons.time_interval import TimeInterval
from commons.utils import addsep
import Sat.SatManager as Sat
from commons.layer import Layer
from commons.mask import Mask
from commons.dataextractor import DataExtractor
from layer_integral.mapbuilder import MapBuilder
import numpy as np
import os
from commons import netcdf3


starttime='20160301'
end__time='20160308'
INPUTDIR='/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data/TMP/'  #args.inputdir
OUTDIR  ='/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data/TMP/'  #args.outdir
maskfile='/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc' # args.maskfile

TI=TimeInterval(starttime,end__time,'%Y%m%d')
archive_dir='/pico/home/usera07ogs/a07ogs00/OPA/V4/archive'
TheMask=Mask(maskfile)

TS = TimeSeries(TI, archive_dir,postfix_dir='POSTPROC/AVE_FREQ_1/',glob_pattern="ave*gz")
forecasts        =TS.get_forecast_days(rundays=[2])
forecasts_sublist=TS.get_sublist(forecasts,[2,3,4]) #forecast tuesday and wed,thu

sat_archive="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/"
DAILY_SAT_LIST=TS.get_daily_sat(forecasts_sublist,sat_archive)

# float aggregator already done by others
day=0
surf_layer=Layer(0,10)
for time,archived_file,satfile in DAILY_SAT_LIST:
    avefile=INPUTDIR + os.path.basename(archived_file)[:-3]
    day=day+1
    outfile=OUTDIR + "misfit+%dh.nc" % (day*24)
    print avefile
    continue
    Sat16   = Sat.convertinV4format( Sat.readfromfile(satfile) )
    De      = DataExtractor(TheMask,filename=avefile, varname='P_i')
    Model   = MapBuilder.get_layer_average(De, surf_layer)

    Misfit = Sat16-Model

    cloudsLand = (np.isnan(Sat16)) #| (Sat16 > 1.e19) | (Sat16<0)
    modelLand  = np.isnan(Model) #lands are nan
    nodata     = cloudsLand | modelLand
    selection  = ~nodata # & TheMask.mask_at_level(200.0)
    Misfit[nodata] = np.NaN
    
    netcdf3.write_2d_file(Misfit, 'chl_misfit', outfile, TheMask)
    
    
    



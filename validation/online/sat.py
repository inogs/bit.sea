from commons.timeseries import TimeSeries
from commons.time_interval import TimeInterval
from commons.utils import addsep

starttime='20160301'
end__time='20160308'
TI=TimeInterval(starttime,end__time,'%Y%m%d')
archive_dir='/pico/home/usera07ogs/a07ogs00/OPA/V4/archive'

TS = TimeSeries(TI, archive_dir,postfix_dir='POSTPROC/AVE_FREQ_1/',glob_pattern="ave*gz")
forecasts=TS.get_forecast_days(rundays=[2])
forecasts_2=TS.get_sublist(forecasts,[2,3]) #forecast tuesday and wed

sat_archive="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/"



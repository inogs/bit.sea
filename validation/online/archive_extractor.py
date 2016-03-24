from commons.timeseries import TimeSeries
from commons.time_interval import TimeInterval

starttime='20160301'
end__time='20160308'
LOC = "/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/online_validation_data/"

TI=TimeInterval(starttime,end__time,'%Y%m%d')

archive_dir="/pico/home/usera07ogs/a07ogs00/OPA/V4/archive/"
T_bio = TimeSeries(TI, archive_dir,postfix_dir='POSTPROC/AVE_FREQ_1/',glob_pattern="ave*gz")
T_phys= TimeSeries(TI, archive_dir,postfix_dir='OPAOPER_A/'          ,glob_pattern="*gz"   )

T_bio.extract_analysis( LOC + 'output_bio/')
T_phys.extract_analysis(LOC + 'output_phys/'); 



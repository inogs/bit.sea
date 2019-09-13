from commons.timeseries import TimeSeries
from commons.time_interval import TimeInterval
from dateutil.relativedelta import relativedelta
TI = TimeInterval("20190801","20190915","%Y%m%d")
archive_dir="/gpfs/work/OGS_prod_0/OPA/V5C/prod/archive/"

LOC="/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/AVE/"
for var in ["P_l", "N3n", "O2o"]:
    T_bio = TimeSeries(TI, archive_dir,postfix_dir='POSTPROC/AVE_FREQ_1/ARCHIVE/',glob_pattern="ave*" +var + ".nc.gz")
    T_bio.extract_analysis(LOC  + "ANALYSISgz/", command = "ln -fs $INFILE $OUTFILE", remove_ext=False)

for var in ["P_l"]:
    T_bio = TimeSeries(TI, archive_dir,postfix_dir='POSTPROC/AVE_FREQ_1/ARCHIVE/',glob_pattern="ave*" +var + ".nc.gz")
    T_bio.extract_forecast(LOC  + "FORECASTgz/", command = "ln -s $INFILE $OUTFILE", remove_ext=False)


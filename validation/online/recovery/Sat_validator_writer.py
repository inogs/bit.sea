from commons.timeseries import TimeSeries
from commons.time_interval import TimeInterval
from dateutil.relativedelta import relativedelta
TI = TimeInterval("20190401","20190801","%Y%m%d")
archive_dir="/gpfs/work/OGS_prod_0/OPA/V5C/prod/archive/"

var="P_l"
LOC="/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/AVE/"
T_bio = TimeSeries(TI, archive_dir,postfix_dir='POSTPROC/AVE_FREQ_1/ARCHIVE/',glob_pattern="ave*" +var + ".nc.gz")



#T_bio.extract_analysis(LOC  + "ANALYSISgz/", command = "ln -s $INFILE $OUTFILE", remove_ext=False)
#T_bio.extract_forecast(LOC  + "FORECASTgz/", command = "ln -s $INFILE $OUTFILE", remove_ext=False)

print "#! /bin/bash"
print "ONLINE_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V5C/"
print "MASKFILE=/gpfs/work/OGS_prod_0/OPA/V5C/prod/wrkdir/2/MODEL/meshmask.nc"
print "SAT_WEEKLY_DIR=${ONLINE_REPO}/SAT/MULTISENSOR/1Km/NRT/WEEKLY_2_24/"
print "SAT_DAILY_DIR=${ONLINE_REPO}/SAT/MULTISENSOR/1Km/NRT/DAILY/CHECKED_24/"
print "FORECAST_DIR=/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/AVE/FORECAST"
print ""

Tuplelist=T_bio.get_runs([2])
for time, path in Tuplelist[1:]:
    RUNDATE=time.strftime("%Y%m%d")
    time_f0=time - relativedelta(days=7)
    time_f1=time - relativedelta(days=6)
    time_f2=time - relativedelta(days=5)
    RUNDATE_DIR="/gpfs/scratch/userexternal/gbolzon0/CHAIN_V5C/archive/" + RUNDATE
    command = "mkdir -p " + RUNDATE_DIR
    print command
    f_0_weekly_name=RUNDATE_DIR + "/Validation_f0_" + time_f0.strftime("%Y%m%d") + "_on_weekly_Sat." + time_f0.strftime("%Y%m%d") + ".nc"
    f_0__daily_name=RUNDATE_DIR + "/Validation_f0_" + time_f0.strftime("%Y%m%d") + "_on_daily_Sat."  + time_f0.strftime("%Y%m%d") + ".nc"
    f_1__daily_name=RUNDATE_DIR + "/Validation_f1_" + time_f0.strftime("%Y%m%d") + "_on_daily_Sat."  + time_f1.strftime("%Y%m%d") + ".nc"
    f_2__daily_name=RUNDATE_DIR + "/Validation_f2_" + time_f0.strftime("%Y%m%d") + "_on_daily_Sat."  + time_f2.strftime("%Y%m%d") + ".nc"
    
    basecommand_w="python SatValidation.py -d %s -f $FORECAST_DIR -s $SAT_WEEKLY_DIR -o %s -m $MASKFILE " 
    basecommand_d="python SatValidation.py -d %s -f $FORECAST_DIR -s $SAT_DAILY_DIR -o %s -m $MASKFILE " 
    
    command_f0_w=basecommand_w  %( time_f0.strftime("%Y%m%d"), f_0_weekly_name)
    command_f0_d=basecommand_d  %( time_f0.strftime("%Y%m%d"), f_0__daily_name)
    command_f1_d=basecommand_d  %( time_f1.strftime("%Y%m%d"), f_1__daily_name)
    command_f2_d=basecommand_d  %( time_f2.strftime("%Y%m%d"), f_2__daily_name)
    for command in [command_f0_w, command_f0_d, command_f1_d, command_f2_d ] :
        print command
    
    
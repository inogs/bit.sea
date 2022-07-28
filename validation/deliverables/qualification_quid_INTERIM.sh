# Descriptor of CMEMS-Med-biogeochemistry-ScQP-...REAN....docx

# QUID INTERIME: 
# SECTION 4:
export MASKFILE=/g100_scratch/userexternal/gcoidess/REANALISI_INTERIM/wrkdir/MODEL/meshmask.nc
# NEW FLOATS DATASET:
#export ONLINE_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V8C/
export ONLINE_REPO=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/DA/new_superfloat/
# INPUT AGGREGATE FREQ2 --> weekly
INPUTDIR=/g100_scratch/userexternal/gcoidess/REANALISI_INTERIM/wrkdir/MODEL/AVE_FREQ_2/
INPUT_AGGR_DIR=/g100_scratch/userexternal/gcoidess/REANALISI_INTERIM/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/
  
# GENERATE THE MAPS:

./plot_MEDEAF_maps.sh


# SAT validation:
mkdir -p Chla_SAT/Tserie Chla_SAT/offshore Chla_SAT/coast Chla_SAT/table/table_CHLA_offshore Chla_SAT/table/table_CHLA_coast 
#SAT_WEEKLY_DIR=/g100_scratch/userexternal/ateruzzi/INTERIM_24_SAT/WEEKLY_4/WEEKLY_4_24/
#SAT_WEEKLY_DIR=/g100_scratch/userexternal/lfeudale/INTERIM_24_SAT/WEEKLY_4/WEEKLY_4_24/
SAT_WEEKLY_DIR=/g100_scratch/userexternal/lfeudale/SAT_INTERIM_24/WEEKLY_4/ALL/WEEKLY_4_24/
#/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/WEEKLY_24_Friday
OUTDIR=TMP_week/
python ScMYvalidation_plan_STD_CORR_valid.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c open_sea   -o export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl
python ScMYvalidation_plan_STD_CORR_valid.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c coast      -o export_data_ScMYValidation_plan_coast_STD_CORR.pkl
#python ScMYvalidation_plan_STD_CORR.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c coast      -o export_data_ScMYValidation_plan_coast_STD_CORR.pkl

##########
python plot_timeseries_STD.py -o export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl -c export_data_ScMYValidation_plan_coast_STD_CORR.pkl -O ./Chla_SAT/Tserie/

# CHL-SURF-W-CLASS4-SIMG-BIAS-BASIN
python plot_timeseries_RMS_CORR.py -i export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl -o Chla_SAT/offshore  # table4.1 ./Chla_SAT/offshore_LONG ./Chla_SAT/coast_LONG
python plot_timeseries_RMS_CORR.py -i export_data_ScMYValidation_plan_coast_STD_CORR.pkl    -o Chla_SAT/coast  

#mv Chla_SAT/offshore/table4.1.dat Chla_Spython plot_timeseries_RMS_CORR.pyAT/table/table_CHLA_offshore.dat
#mv Chla_SAT/coast/table4.1.dat Chla_SAT/table/table_CHLA_coast.dat

import sys
sys.exit()
# BIOFLOATS SECTION: Hovmoeller plots, wmo trajectories and statistics per basin
BASEDIR=/g100_scratch/userexternal/lfeudale/REANALYSIS/REAN24/validation/INTERIM2021/PROFILATORE/

NCDIR=tmp_nc
OUTDIR=Hovmoeller_plots #Fig4.4

mkdir $NCDIR $OUTDIR
echo python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -b $BASEDIR -o $NCDIR
python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -b $BASEDIR -o $NCDIR
echo python Hov_Stat_plot.py -m $MASKFILE -i $NCDIR -o $OUTDIR -b $BASEDIR
python Hov_Stat_plot.py -m $MASKFILE -i $NCDIR -o $OUTDIR -b $BASEDIR

#OUTDIR=Fig4.6
#OUTDIR=FLOATVAR-PROF-D-CLASS4-PROF-CORR-BASIN
##mkdir -p $NCDIR/chla_check/ $NCDIR/nitrate_check/
#python BASIN_Float_vs_Model_Stat_Timeseries_monthly.py -m $MASKFILE -o $NCDIR
#python BASIN_Float_vs_Model_Stat_Timeseries_monthly_plotter.py -m $MASKFILE -i $NCDIR -o $OUTDIR
#mv $OUTDIR/N3n*.png Fig4.13


OUTFIGDIR=Floats_bias_rmse_Timeseries     # 8layer x 7sub x 3var = 168 png files
TABLE_DIR=Floats_bias_rmse_tables         #: 2stats x 3var        = 6 txt files, TABLE.O2o_BIAS.txt  with time average for each layer,sub
mkdir -p $OUTFIGDIR $TABLE_DIR #table4.3/ table4.9/ table4.12/ Fig4.5/ Fig4.14/ Fig4.16/
python biofloats_ms.py  -m $MASKFILE -o float_bias_rmse.nc
python biofloats_ms_plotter.py -i float_bias_rmse.nc -f $OUTFIGDIR -t $TABLE_DIR


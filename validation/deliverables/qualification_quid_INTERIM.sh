# Descriptor of CMEMS-Med-biogeochemistry-ScQP-...REAN....docx

# QUID INTERIME: 
# SECTION 4:
export MASKFILE=/g100_scratch/userexternal/gcoidess/NEW_REA_24/wrkdir/MODEL/meshmask.nc

# NEW FLOATS DATASET:
export ONLINE_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V9C/
# INPUT AGGREGATE FREQ2 --> weekly
INPUTDIR=/g100_scratch/userexternal/gcoidess/NEW_REA_24/wrkdir/MODEL/AVE_FREQ_2/
INPUT_AGGR_DIR=/g100_scratch/userexternal/gcoidess/NEW_REA_24/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/
  
# GENERATE THE MAPS:

./plot_MEDEAF_maps.sh


# SAT validation:
mkdir -p Chla_SAT/Tserie Chla_SAT/offshore Chla_SAT/coast Chla_SAT/table/table_CHLA_offshore Chla_SAT/table/table_CHLA_coast 
SAT_WEEKLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V9C/SAT/CHL/DT/WEEKLY_4_24/

OUTDIR=tmp_WEEK/
mkdir $OUTDIR
echo python ScMYvalidation_plan.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c open_sea   -o export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl -v "chl" -l 10
python ScMYvalidation_plan_STD_CORR_valid.py  -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c open_sea   -o export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl # -v "chl" -l 10
python ScMYvalidation_plan_STD_CORR_valid.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c coast      -o export_data_ScMYValidation_plan_coast_STD_CORR.pkl # -v "chl" -l 10

##########
python plot_timeseries_STD.py -v chl -i export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl -o ./Chla_SAT/Tserie/ 

# CHL-SURF-W-CLASS4-SIMG-BIAS-BASIN
python plot_timeseries_RMS_CORR.py -v chl -i export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl -o Chla_SAT/offshore  # table4.1 ./Chla_SAT/offshore_LONG ./Chla_SAT/coast_LONG
python plot_timeseries_RMS_CORR.py -v chl -i export_data_ScMYValidation_plan_coast_STD_CORR.pkl    -o Chla_SAT/coast  


exit 0

# BIOFLOATS SECTION: Hovmoeller plots, wmo trajectories and statistics per basin
BASEDIR=/g100_scratch/userexternal/lfeudale/REANALYSIS/REAN24/validation/NEW_REA_2020_2022/PROFILATORE/

NCDIR=tmp_nc
OUTDIR=Hovmoeller_plots #Fig4.4

mkdir $NCDIR $OUTDIR
echo python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -b $BASEDIR -o $NCDIR
python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -b $BASEDIR -o $NCDIR
echo python Hov_Stat_plot.py -m $MASKFILE -i $NCDIR -o $OUTDIR -b $BASEDIR
python Hov_Stat_plot.py -m $MASKFILE -i $NCDIR -o $OUTDIR -b $BASEDIR


OUTFIGDIR=Floats_bias_rmse_Timeseries     # 8layer x 7sub x 3var = 168 png files
TABLE_DIR=Floats_bias_rmse_tables         #: 2stats x 3var        = 6 txt files, TABLE.O2o_BIAS.txt  with time average for each layer,sub
mkdir -p $OUTFIGDIR $TABLE_DIR #table4.3/ table4.9/ table4.12/ Fig4.5/ Fig4.14/ Fig4.16/
python biofloats_ms.py  -m $MASKFILE -o float_bias_rmse.nc
python biofloats_ms_plotter.py -i float_bias_rmse.nc -f $OUTFIGDIR -t $TABLE_DIR


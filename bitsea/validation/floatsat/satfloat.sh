# from qualification_quid.sh in ../deliverables

export RUN=DA_FLOAT_SAT/Summer/RUN_FLOAT_02/
export MASKFILE=/pico/scratch/userexternal/ateruzzi/$RUN/wrkdir/MODEL/meshmask.nc


mkdir -p $RUN/Fig4.2 $RUN/Fig4.3/offshore $RUN/table4.1 $RUN/table4.2
SAT_DAILY_DIR=/pico/scratch/userexternal/ateruzzi/SATELLITE/MULTISENSOR/DAILY/CHECKED_INTERP/
INPUT_AGGR_DIR=/pico/scratch/userexternal/ateruzzi/$RUN/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/

python ScMYvalidation_plan.py -s $SAT_DAILY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c open_sea   -o $RUN/export_data_ScMYValidation_plan_open_sea.pkl
python plot_timeseries_2014_2016_openonly.py -m $MASKFILE -o $RUN/export_data_ScMYValidation_plan_open_sea.pkl -O $RUN/Fig4.2/

python plot_timeseries_RMS_statsall.py -i $RUN/export_data_ScMYValidation_plan_open_sea.pkl -o $RUN/Fig4.3/offshore  # table4.1
cp $RUN/Fig4.3/offshore/table4.1.dat $RUN/table4.1


# BIOFLOATS SECTION: Hovmoeller plots, wmo trajectories and statistics per basin
# The first time execute profiler_floatsat.py
BASEDIR=$RUN/PROFILATORE/
mkdir -p $BASEDIR
python profiler_floatsat.py #to be executed one time
# ATTENTION: profiler is imported in the following scripts.
#            Use consistent run name


# CHECK THE RUN NAME IN profiler_floatsat.py!
# Figures 4.4a
mkdir -p $RUN/Fig4.4a $RUN/Fig4.4b $RUN/Fig4.6 $RUN/tmp_nc $RUN/table4.4
mkdir -p $RUN/PROFILATORE
export OUTDIR=$RUN/Fig4.4a
#python Hov_flots+model.py -m $MASKFILE -o $OUTDIR
# MODIFY the dates in the script Hov_flots+model_vars_runs.py
python Hov_flots+model_vars_runs.py -m $MASKFILE -o $OUTDIR -r $RUN
#mv $OUTDIR/*N3n*.png Fig4.12
#mv $OUTDIR/*O2o*.png Fig4.15


# Figures 4.4b
NCDIR=$RUN/tmp_nc
OUTDIR=$RUN/Fig4.4b
mkdir -p $NCDIR
mkdir -p $OUTDIR
python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -o $NCDIR
python SingleFloat_vs_Model_Stat_Timeseries_plotter.py -i $NCDIR -o $OUTDIR
#mv $OUTDIR/N3n*.png Fig4.12
#mv $OUTDIR/O2o*.png Fig4.15

# Figures 4.6 and tables 4.4 
# subbasin DCM and MLB
NCDIR=$RUN/tmp_nc
OUTDIR=$RUN/Fig4.6
mkdir -p $OUTDIR $RUN/table4.4
python BASIN_Float_vs_Model_Stat_Timeseries_monthly.py -m $MASKFILE -o $NCDIR
python BASIN_Float_vs_Model_Stat_Timeseries_monthly_plotter.py -m $MASKFILE -i $NCDIR -o $OUTDIR
#mv $OUTDIR/N3n*.png Fig4.13

cp $OUTDIR/P_l_tab_statistics_SHORT.txt $RUN/table4.4/ 
#cp $OUTDIR/N3n_tab_statistics_SHORT.txt table4.8/


# Comparison on DA dates (not on float dates)
# The first time execute profiler_floatsat_DAdates.py
BASEDIR=$RUN/PROFILATORE_DAtimes/
mkdir -p $BASEDIR
python profiler_floatsat_DAdates.py #to be executed one time
# ATTENTION: profiler is imported in the following scripts.
#            Use consistent run name


# Statistics on DA dates
NCDIR=$RUN/tmp_nc_DAdates
mkdir -p $NCDIR 
python SingleFloat_vs_Model_Stat_TS_DAdates.py -m $MASKFILE -o $NCDIR
# Plots with all runs in allruns.sh



exit 0

# BIOFLOATS SECTION: statistics on layers
# CHL-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN
# NIT-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN
#  DO-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN


OUTFIGDIR=Floats_bias_rmse_Timeseries     # 8layer x 7sub x 3var = 168 png files
TABLE_DIR=Floats_bias_rmse_tables         #: 2stats x 3var        = 6 txt files, TABLE.O2o_BIAS.txt  with time average for each layer,sub
mkdir -p $OUTFIGDIR $TABLE_DIR table4.3/ table4.9/ table4.12/ Fig4.5/ Fig4.14/ Fig4.16/
python biofloats_ms.py  -m $MASKFILE -o float_bias_rmse.nc
python biofloats_ms_plotter.py -i float_bias_rmse.nc -f $OUTFIGDIR -t $TABLE_DIR
cp $TABLE_DIR/P_l_BIAS.txt $TABLE_DIR/P_l_RMSE.txt table4.3/
cp $TABLE_DIR/N3n_BIAS.txt $TABLE_DIR/N3n_RMSE.txt table4.9/
cp $TABLE_DIR/O2o_BIAS.txt $TABLE_DIR/O2o_RMSE.txt table4.12/
cp $OUTFIGDIR/*P_l* Fig4.5
cp $OUTFIGDIR/*N3n* Fig4.14
cp $OUTFIGDIR/*O2o* Fig4.16
#########################   static dataset climatology section ###################################
# Fig. 4.11:
# PHO-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers
# NIT-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers
#  DO-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers

# Fig. 4.19:
# ALK-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers
# DIC-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers

# Figures 4.11 and 4.19
mkdir -p sim_vs_clim_profiles/ Fig4.11 Fig4.19
python simulation_vs_clim.py -i $STAT_PROFILES_DIR -o sim_vs_clim_profiles/ -s 20150101 -e 20170101 -m $MASKFILE
cp sim_vs_clim_profiles/Fig_4.11*png Fig4.11
cp sim_vs_clim_profiles/Fig_4.19*png Fig4.19

DIR=static_clim
mkdir -p $DIR table4.6 table4.7 table4.9/ table4.11 table4.13/ table4.14
# -------------------------------------------------------------------------

python static_clim_validation.py -i $STAT_PROFILES_DIR -o $DIR -m $MASKFILE -s 20150101 -e 20170101

# -------------------------------------------------------------------------
cp $DIR/N1p-LAYER-Y-CLASS4-CLIM.txt $DIR/N3n-LAYER-Y-CLASS4-CLIM.txt table4.6/
cp $DIR/O2o-LAYER-Y-CLASS4-CLIM.txt                                  table4.10/
cp $DIR/Ac-LAYER-Y-CLASS4-CLIM.txt $DIR/DIC-LAYER-Y-CLASS4-CLIM.txt  table4.13/

# ALK-PROF-Y-CLASS4-CLIM-CORR-BASIN calculated on 14 layers --> in table4.14
# DIC-PROF-Y-CLASS4-CLIM-CORR-BASIN calculated on 14 layers --> in table4.14

cp $DIR/N1p-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt $DIR/N3n-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt table4.7
cp $DIR/O2o-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt                                            table4.11
cp $DIR/Ac-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt  $DIR/DIC-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt table4.14

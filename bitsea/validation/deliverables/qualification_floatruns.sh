
export MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS16corrected/meshmask.nc


RUN=DA_Float/RUN_REF/
RUN=DA_Float/RUN_SAT12/
RUN=DA_Float/RUN_FLOAT_chl12/
RUN=DA_Float/RUN_FLOAT_chl_n/
RUN=DA_Float/RUN_FLOAT_chl_nupd/


   INPUT_AGGR_DIR=$CINECA_SCRATCH/$RUN/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/
SAT_WEEKLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MULTISENSOR/1Km/DT/WEEKLY_2_16/
OUTDIR=/gpfs/scratch/userexternal/ateruzzi/ELAB_DAFloat/QuIDeval/$RUN/
mkdir -p $OUTDIR
echo python ScMYvalidation_plan_STD_CORR.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c open_sea   -o $OUTDIR/export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl
# python ScMYvalidation_plan_STD_CORR.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c coast      -o $OUTDIR/export_data_ScMYValidation_plan_coast_STD_CORR.pkl


#python plot_timeseries_STD.py -o export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl -c export_data_ScMYValidation_plan_coast_STD_CORR.pkl -O ./Fig4.2/

# CHL-SURF-W-CLASS4-SIMG-BIAS-BASIN
# CHL-SURF-W-CLASS4-SIMG-RMSD-BASIN
OUTFIG=$OUTDIR/Fig4.3/
mkdir -p $OUTFIG/offshore
echo python plot_timeseries_RMS_CORR.py -i $OUTDIR/export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl -o $OUTFIG/offshore  # table4.1
#python plot_timeseries_RMS_CORR.py -i export_data_ScMYValidation_plan_coast_STD_CORR.pkl    -o Fig4.3/coast     # table4.2
#cp Fig4.3/offshore/table4.1.dat table4.1
#cp Fig4.3/coast/table4.1.dat    table4.2/table4.2.dat


# Execute profiler.py once for each run
export ONLINE_REPO=/gpfs/scratch/userexternal/ateruzzi/REPO_LOV/ONLINE/ # THE SAME OF DA

# BIOFLOATS SECTION: Hovmoeller plots, wmo trajectories and statistics per basin
# It is included the 2014 also.
# Figures 4.4a
mkdir -p $OUTDIR/Fig4.4a $OUTDIR/Fig4.4b $OUTDIR/Fig4.6 $OUTDIR/Fig4.12 $OUTDIR//Fig4.13 $OUTDIR/Fig4.15 $OUTDIR/tmp_nc $OUTDIR/table4.4 $OUTDIR/table4.8 
OUTFIG=$OUTDIR/Fig4.4a
echo python Hov_flots+model_vars.py -m $MASKFILE -o $OUTFIG
mv $OUTFIG/*N3n*.png $OUTDIR/Fig4.12
mv $OUTFIG/*O2o*.png $OUTDIR/Fig4.15

# Figures 4.4b
NCDIR=$OUTDIR/tmp_nc
OUTFIG=$OUTDIR/Fig4.4b
echo python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -o $NCDIR
echo python SingleFloat_vs_Model_Stat_Timeseries_plotter.py -i $NCDIR -o $OUTFIG -m $MASKFILE
mv $OUTFIG/N3n*.png $OUTDIR/Fig4.12
mv $OUTFIG/O2o*.png $OUTDIR/Fig4.15

# Figures 4.6 and 4.13 + tables 4.4 and 4.8
# CHL-PROF-D-CLASS4-PROF-CORR-BASIN
# NIT-PROF-D-CLASS4-PROF-CORR-BASIN
#  DO-PROF-D-CLASS4-PROF-CORR-BASIN 

OUTFIG=$OUTDIR/Fig4.6
echo python BASIN_Float_vs_Model_Stat_Timeseries_monthly.py -m $MASKFILE -o $NCDIR
echo python BASIN_Float_vs_Model_Stat_Timeseries_monthly_plotter.py -m $MASKFILE -i $NCDIR -o $OUTFIG
mv $OUTFIG/N3n*.png $OUTDIR/Fig4.13

cp $OUTFIG/P_l_tab_statistics_SHORT.txt $OUTDIR/table4.4/ 
cp $OUTFIG/N3n_tab_statistics_SHORT.txt $OUTDIR/table4.8/



# BIOFLOATS SECTION: statistics on layers
# CHL-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN
# NIT-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN
#  DO-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN


OUTFIGDIR=$OUTDIR/Floats_bias_rmse_Timeseries     # 8layer x 7sub x 3var = 168 png files
TABLE_DIR=$OUTDIR/Floats_bias_rmse_tables         #: 2stats x 3var        = 6 txt files, TABLE.O2o_BIAS.txt  with time average for each layer,sub
mkdir -p $OUTFIGDIR $TABLE_DIR $OUTDIR/table4.3/ $OUTDIR/table4.9/ $OUTDIR/table4.12/ $OUTDIR/Fig4.5/ $OUTDIR/Fig4.14/ $OUTDIR/Fig4.16/
echo python biofloats_ms.py  -m $MASKFILE -o $OUTDIR/float_bias_rmse.nc
echo python biofloats_ms_plotter.py -i $OUTDIR/float_bias_rmse.nc -f $OUTFIGDIR -t $TABLE_DIR
cp $TABLE_DIR/P_l_BIAS.txt $TABLE_DIR/P_l_RMSE.txt $OUTDIR/table4.3/
cp $TABLE_DIR/N3n_BIAS.txt $TABLE_DIR/N3n_RMSE.txt $OUTDIR/table4.9/
cp $TABLE_DIR/O2o_BIAS.txt $TABLE_DIR/O2o_RMSE.txt $OUTDIR/table4.12/
cp $OUTFIGDIR/*P_l* $OUTDIR/Fig4.5
cp $OUTFIGDIR/*N3n* $OUTDIR/Fig4.14
cp $OUTFIGDIR/*O2o* $OUTDIR/Fig4.16
#########################   static dataset climatology section ###################################
# Fig. 4.11:
# PHO-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers
RUN=DA_Float/RUN_FLOAT_chl_n/
n biofloats_ms_plotter.py -i /gpfs/scratch/userexternal/ateruzzi/ELAB_DAFloat/QuIDeval/DA_Float/RUN_FLOAT_chl_n///float_bias_rmse.nc -f /gpfs/scratch/userexternal/ateruzzi/ELAB_DAFloat/QuIDeval/DA_Float/RUN_FLOAT_chl_n///Floats_bias_rmse_Timeseries -t /gpfs/scratch/userexternal/ateruzzi/ELAB_DAFloat/QuIDeval/DA_Float/RUN_FLOAT_chl_n///Floats_bias_rmse_tables

# NIT-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers
#  DO-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers

# Fig. 4.19:
# ALK-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers
# DIC-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers
 
# Fig. 4.17:
# PH-PROF-Y-CLASS4-[CLIM/LIT]-MEAN
# pCO-PROF-Y-CLASS4-[CLIM/LIT]-MEAN

#STAT_PROFILES_MONTHLY_DIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v20_1/wrkdir/POSTPROC/output/MONTHLY/STAT_PROFILES/
STAT_PROFILES_MONTHLY_DIR=/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_05/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/
# Figures 4.11 and 4.17 (previously named 4.19)
mkdir -p sim_vs_clim_profiles/ Fig4.11 Fig4.17
python simulation_vs_clim.py -i $STAT_PROFILES_MONTHLY_DIR -o sim_vs_clim_profiles/ -s 20170101 -e 20180101 -m $MASKFILE
cp sim_vs_clim_profiles/Fig_4.11*png Fig4.11
cp sim_vs_clim_profiles/Fig_4.17*png Fig4.17

DIR=static_clim
mkdir -p $DIR table4.6 table4.7 table4.10/ table4.11 table4.13/ table4.14
# -------------------------------------------------------------------------

#python static_clim_validation_STD.py -i $STAT_PROFILES_DIR -o $DIR -m $MASKFILE -s 20170101 -e 20180101
#python static_clim_validation_STD_pH.py -i $STAT_PROFILES_DIR -o $DIR -m $MASKFILE -s 20170101 -e 20180101
python static_clim_validation.py -i $STAT_PROFILES_DIR -o $DIR -m $MASKFILE -s 20170101 -e 20180101

# -------------------------------------------------------------------------
cp $DIR/N1p-LAYER-Y-CLASS4-CLIM.txt $DIR/N3n-LAYER-Y-CLASS4-CLIM.txt table4.6/
cp $DIR/O2o-LAYER-Y-CLASS4-CLIM.txt                                  table4.10/

# PH-PROF-Y-CLASS4-CLIM-[BIAS/RMSD/CORR]-BASIN
# PCO-PROF-Y-CLASS4-CLIM-[BIAS/RMSD/CORR]-BASIN

#cp $DIR/Ac-LAYER-Y-CLASS4-CLIM.txt $DIR/DIC-LAYER-Y-CLASS4-CLIM.txt  table4.13/
cp $DIR/ALK-LAYER-Y-CLASS4-CLIM.txt $DIR/DIC-LAYER-Y-CLASS4-CLIM.txt  table4.13/
cp $DIR/pH-LAYER-Y-CLASS4-CLIM.txt  $DIR/pCO2-LAYER-Y-CLASS4-CLIM.txt table4.13/

# ALK-PROF-Y-CLASS4-CLIM-CORR-BASIN calculated on 14 layers --> in table4.14
# DIC-PROF-Y-CLASS4-CLIM-CORR-BASIN calculated on 14 layers --> in table4.14

cp $DIR/N1p-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt $DIR/N3n-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt table4.7
cp $DIR/O2o-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt                                            table4.11
#cp $DIR/Ac-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt  $DIR/DIC-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt table4.14
cp $DIR/ALK-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt  $DIR/DIC-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt table4.14
cp $DIR/pH-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt   $DIR/pCO2-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt table4.14

### SECTION WITH NEW FIGURES FOR pCO2 COMPARISON/ANALYSIS ###
# Compare pCO2 vs SOCAT dataset: 
# Fig4.20  : PCO-SURF-M-CLASS4-CLIM-MEAN-BASIN
# table4.15: PCO-SURF-M-CLASS4-CLIM-RMSD-BASIN.txt 
mkdir -p table4.15 Fig4.20
# First generate climatology tables:
python monthly_clim_socat_pCO2.py
#python monthly_2017.py

mkdir -p monthly_2017_surf/
python monthly_2017 -i $STAT_PROFILES_DIR -o monthly_2017_surf/

python table_pCO2vsSOCAT.py -i monthly_2017_surf/
mv pCO2-SURF-M-CLASS4-CLIM-RMSD-BASIN.txt table4.15/.
mv TOT_RMSD_pCO2vsSOCAT.txt table4.15/.

# generate new figure 4.20: comparison pCO2 BFM vs SOCAT dataset

python plot_month_pCO2vsSOCAT.py
mv pCO2_monthly_tseries_Fig4.20.png Fig4.20/.

# generate new figure 4.21: comparison pCO2 and CO2ariflux vs REAN (CLIM)
# Fig4.21: airseaCO2flux-SURF-M-CLASS4-CLIM-MEAN-BASIN

REANDIR="/gpfs/scratch/userexternal/gbolzon0/REA/output/STAT_PROFILES/"
#python plot_pCO2_CO2airflux_vs_REAN.py -r "/gpfs/scratch/userexternal/gbolzon0/REA/output/STAT_PROFILES/"
python plot_pCO2_CO2airflux_vs_REAN.py -r $REANDIR -i $STAT_PROFILES_MONTHLY_DIR -o Fig4.21/


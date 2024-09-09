# Descriptor of CMEMS-Med-biogeochemistry-ScQP-v1.3.docx

# QUID REANALYSIS
# SECTION 4:
export MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MULTIVARIATE_24/TEST_02/wrkdir/MODEL/meshmask.nc
export ONLINE_REPO=/gpfs/scratch/userexternal/ateruzzi/REPO_ALL_FLOATS/ONLINE_V5C/

         INPUTDIR=/gpfs/scratch/userexternal/ateruzzi/MULTIVARIATE_24/TEST_02/wrkdir/MODEL/AVE_FREQ_2/ #compressed/
   INPUT_AGGR_DIR=/gpfs/scratch/userexternal/ateruzzi/MULTIVARIATE_24/TEST_02/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/
STAT_PROFILES_DIR=/gpfs/scratch/userexternal/ateruzzi/MULTIVARIATE_24/TEST_02/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/
SAT_MONTHLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/MONTHLY_24/

OUTDIR=spaghettiplots
mkdir -p $OUTDIR
python plot_layer_timeseries_on_profiles.py -i $STAT_PROFILES_DIR -m $MASKFILE -o $OUTDIR

mkdir -p Fig4.1/ Fig4.7 Fig4.9 Fig4.10 Fig4.18 Fig4.19 Fig4.7bis


COMMONS_PARAMS="-m $MASKFILE -o LayerMaps/  -l Plotlist_bio.xml -s 20170101 -e 20180101"

mkdir -p LayerMaps
python averager_and_plot_map.py -i $INPUT_AGGR_DIR  -v P_l  -t mean $COMMONS_PARAMS      # Fig4.1 CHL-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v N3n  -t mean $COMMONS_PARAMS      # Fig4.10 NIT-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v N1p  -t mean $COMMONS_PARAMS      # Fig4.9  PHOS-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
#INPUTDIR_PPN=/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_05/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/
#python averager_and_plot_map_ppn.py -i $INPUTDIR_PPN        -v ppn  -t integral $COMMONS_PARAMS  # Fig4.7 NPP-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN per lo 0-200m
#python averager_and_plot_map_ppn.py -i $INPUTDIR_PPN        -v ppn  -t integral -m $MASKFILE -o Fig4.7bis -l Plotlist_bio_ppn20m.xml -s 20170101 -e 20180101

#INPUTDIR_PPN=/gpfs/scratch/userexternal/ateruzzi/PPN_TRANSITION2017/
mkdir -p Fig4.7refScale
#python averager_and_plot_map_ppn_refScale.py -i $INPUTDIR_PPN  -v ppn  -t integral -m $MASKFILE -o Fig4.7refScale -l Plotlist_bio.xml -s 20170101 -e 20180101
python averager_and_plot_map_ppn_refScale.py -i $INPUT_AGGR_DIR  -v ppn  -t integral -m $MASKFILE -o Fig4.7refScale -l Plotlist_bio.xml -s 20170101 -e 20180101

python averager_and_plot_map.py -i $INPUTDIR        -v ALK   -t mean $COMMONS_PARAMS   # Ac-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN  --> not requested in the ScQP
python averager_and_plot_map.py -i $INPUTDIR        -v DIC  -t mean $COMMONS_PARAMS   # DIC-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN --> not requested in the ScQP
cp LayerMaps/*P_l* Fig4.1/
cp LayerMaps/*N3n* Fig4.10/
cp LayerMaps/*N1p* Fig4.9/
cp LayerMaps/*ppn* Fig4.7/
cp LayerMaps/*ALK* Fig4.18/
cp LayerMaps/*DIC* Fig4.19/

# python averager_and_plot_map.py -i $INPUTDIR        -v O2o  -t mean $COMMONS_PARAMS  # DO-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN --> Not showed for lack of ref
# python averager_and_plot_map.py -i $INPUT_AGGR_DIR  -v P_c  -t mean $COMMONS_PARAMS  # PHYC-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN --> Not showed for lack of ref

#CHL-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN from SATELLITE:
python sat_ave_and_plot.py      -i $SAT_MONTHLY_DIR -m $MASKFILE  -o Fig4.1/

mkdir -p Fig4.2 Fig4.3/offshore Fig4.3/coast table4.1 table4.2
SAT_WEEKLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/WEEKLY_24_Friday
###STAT_PROFILES_DIR1=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v19_3/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/
OUTDIR=TMP_week/
#python ScMYvalidation_plan.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c open_sea   -o export_data_ScMYValidation_plan_open_sea.pkl
#python ScMYvalidation_plan.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c coast      -o export_data_ScMYValidation_plan_coast.pkl
python ScMYvalidation_plan_STD_CORR.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c open_sea   -o export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl
python ScMYvalidation_plan_STD_CORR.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c coast      -o export_data_ScMYValidation_plan_coast_STD_CORR.pkl

#python plot_timeseries_2014_2016.py -i $STAT_PROFILES_DIR1 -m $MASKFILE -o export_data_ScMYValidation_plan_open_sea.pkl -c export_data_ScMYValidation_plan_coast.pkl -O ./Fig4.2/
python plot_timeseries_STD.py -o export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl -c export_data_ScMYValidation_plan_coast_STD_CORR.pkl -O ./Fig4.2/

# CHL-SURF-W-CLASS4-SIMG-BIAS-BASIN
# CHL-SURF-W-CLASS4-SIMG-RMSD-BASIN
python plot_timeseries_RMS_CORR.py -i export_data_ScMYValidation_plan_open_sea_STD_CORR.pkl -o Fig4.3/offshore  # table4.1
python plot_timeseries_RMS_CORR.py -i export_data_ScMYValidation_plan_coast_STD_CORR.pkl    -o Fig4.3/coast     # table4.2
cp Fig4.3/offshore/table4.1.dat table4.1
cp Fig4.3/coast/table4.1.dat    table4.2/table4.2.dat

mkdir -p Fig4.8
INTEGRALS_PPN=/gpfs/scratch/userexternal/ateruzzi/MULTIVARIATE_24/TEST_02/wrkdir/POSTPROC/output/PPN_MONTHLY/INTEGRALS/PPN/
python read_ppn_from_avescan_do_plot.py -c open_sea   -i $INTEGRALS_PPN -o Fig4.8


# Figure carbonatiche 1x1
# ALK-LAYER-Y-CLASS4-CLIM-RMSD
# ALK-LAYER-Y-CLASS4-CLIM-BIAS
# DIC-LAYER-Y-CLASS4-CLIM-RMS
# DIC-LAYER-Y-CLASS4-CLIM-BIAS

# ---------- executed elsewhere
# python ricostruzione_Integrals.py -i /pico/scratch/userexternal/gbolzon0/eas_v12/eas_v19_3/wrkdir/POSTPROC/output/AVE_FREQ_2/1x1/INTEGRALS/ -o 1x1/
# ----------------------------
DIR_MOD_1x1=/galileo/home/userexternal/lfeudale/matlab_MAP_1x1/2019/
DIR_CLI_1x1=/galileo/home/userexternal/lfeudale/matlab_MAP_1x1/MAP1x1/
DIR_OUT=.
python readMAP1x1_8layer_do_CFR_carbsys.py -i $DIR_MOD_1x1 -c $DIR_CLI_1x1 -o $DIR_OUT

#### PH-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN in Fig4.20
# OLD FIGURE:
mkdir -p Fig4.20/Feb Fig4.20/May Fig4.20/Aug Fig4.20/Nov Fig4.21/Feb Fig4.21/May  Fig4.21/Aug  Fig4.21/Nov
INPUTDIR_MON=/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_05/wrkdir/MODEL/AVE_FREQ_2/
COMMONS_PARAMS="-i $INPUTDIR_MON -m $MASKFILE -l Plotlist_bio.xml  -t mean "
python averager_and_plot_map.py $COMMONS_PARAMS  -v pH  -s 20170201 -e 20170301 -o Fig4.20/Feb
python averager_and_plot_map.py $COMMONS_PARAMS  -v pH  -s 20170501 -e 20170601 -o Fig4.20/May
python averager_and_plot_map.py $COMMONS_PARAMS  -v pH  -s 20170801 -e 20170901 -o Fig4.20/Aug
python averager_and_plot_map.py $COMMONS_PARAMS  -v pH  -s 20171101 -e 20171201 -o Fig4.20/Nov

# PCO-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN in Fig4.21
# OLD FIGURE:
python averager_and_plot_map.py $COMMONS_PARAMS  -v pCO2  -s 20170201 -e 20170301 -o Fig4.21/Feb
python averager_and_plot_map.py $COMMONS_PARAMS  -v pCO2  -s 20170501 -e 20170601 -o Fig4.21/May
python averager_and_plot_map.py $COMMONS_PARAMS  -v pCO2  -s 20170801 -e 20170901 -o Fig4.21/Aug
python averager_and_plot_map.py $COMMONS_PARAMS  -v pCO2  -s 20171101 -e 20171201 -o Fig4.21/Nov

# BIOFLOATS SECTION: Hovmoeller plots, wmo trajectories and statistics per basin
# It is included the 2014 also.
# Figures 4.4a
#mkdir -p Fig4.4a Fig4.4b Fig4.6 Fig4.12 Fig4.13 Fig4.15 tmp_nc table4.4 table4.8 
mkdir -p Fig4.4 Fig4.6 Fig4.12 Fig4.13 Fig4.15 tmp_nc table4.4 table4.8
BASEDIR=/gpfs/scratch/userexternal/lfeudale/validation/TEST_02/PROFILATORE/
NCDIR=tmp_nc
OUTDIR=Fig4.4
echo python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -b $BASEDIR -o $NCDIR
python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -b $BASEDIR -o $NCDIR
echo python Hov_Stat_plot.py -m $MASKFILE -i $NCDIR -o $OUTDIR -b $BASEDIR
python Hov_Stat_plot.py -m $MASKFILE -i $NCDIR -o $OUTDIR -b $BASEDIR
#python Hov_flots+model_vars.py -m $MASKFILE -o $OUTDIR
mv $OUTDIR/*N3n*.png Fig4.12
mv $OUTDIR/*O2o*.png Fig4.15

# Figures 4.4b
#NCDIR=tmp_nc
#OUTDIR=Fig4.4b
#python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -o $NCDIR
#python SingleFloat_vs_Model_Stat_Timeseries_plotter.py -i $NCDIR -o $OUTDIR -m $MASKFILE
#mv $OUTDIR/N3n*.png Fig4.12
#mv $OUTDIR/O2o*.png Fig4.15

# Figures 4.6 and 4.13 + tables 4.4 and 4.8
# CHL-PROF-D-CLASS4-PROF-CORR-BASIN
# NIT-PROF-D-CLASS4-PROF-CORR-BASIN
#  DO-PROF-D-CLASS4-PROF-CORR-BASIN 

OUTDIR=Fig4.6
#mkdir -p $NCDIR/chla_check/ $NCDIR/nitrate_check/
python BASIN_Float_vs_Model_Stat_Timeseries_monthly.py -m $MASKFILE -o $NCDIR
python BASIN_Float_vs_Model_Stat_Timeseries_monthly_plotter.py -m $MASKFILE -i $NCDIR -o $OUTDIR
mv $OUTDIR/N3n*.png Fig4.13

cp $OUTDIR/P_l_tab_statistics_SHORT.txt table4.4/ 
cp $OUTDIR/N3n_tab_statistics_SHORT.txt table4.8/



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
 
# Fig. 4.17:
# PH-PROF-Y-CLASS4-[CLIM/LIT]-MEAN
# pCO-PROF-Y-CLASS4-[CLIM/LIT]-MEAN

#STAT_PROFILES_MONTHLY_DIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v20_1/wrkdir/POSTPROC/output/MONTHLY/STAT_PROFILES/
STAT_PROFILES_MONTHLY_DIR=/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_05/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/
# Figures 4.11 and 4.17 (previously named 4.19)
mkdir -p sim_vs_clim_profiles/ Fig4.11 Fig4.17
#python simulation_vs_clim.py -i $STAT_PROFILES_MONTHLY_DIR -o sim_vs_clim_profiles/ -s 20170101 -e 20180101 -m $MASKFILE
python simulation_vs_clim.py -i $STAT_PROFILES_DIR -o sim_vs_clim_profiles/ -s 20170101 -e 20180101 -m $MASKFILE
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
mkdir -p table4.15 Fig4.20 Fig4.21
# First generate climatology tables:
python monthly_clim_socat_pCO2.py
#python monthly_2017.py

mkdir -p monthly_2017_surf/
python monthly_2017.py -i $STAT_PROFILES_DIR -o monthly_2017_surf/

python table_pCO2vsSOCAT.py -i monthly_2017_surf/
mv pCO2-SURF-M-CLASS4-CLIM-RMSD-BASIN.txt table4.15/.
mv TOT_RMSD_pCO2vsSOCAT.txt table4.15/.

# generate new figure 4.20: comparison pCO2 BFM vs SOCAT dataset

python plot_month_pCO2vsSOCAT.py
mv pCO2_monthly_tseries_Fig4.20.png Fig4.20/.

# generate new figure 4.21: comparison pCO2 and CO2ariflux vs REAN (CLIM)
# Fig4.21: airseaCO2flux-SURF-M-CLASS4-CLIM-MEAN-BASIN

#REANDIR="/gpfs/scratch/userexternal/gbolzon0/REA/output/STAT_PROFILES/"
REANDIR=/gpfs/scratch/userexternal/lfeudale/REANALYSIS/STAT_PROFILES/
#python plot_pCO2_CO2airflux_vs_REAN.py -r "/gpfs/scratch/userexternal/gbolzon0/REA/output/STAT_PROFILES/"
python plot_pCO2_CO2airflux_vs_REAN.py -r $REANDIR -i $STAT_PROFILES_MONTHLY_DIR -o Fig4.21/
python plot_pCO2_CO2airflux_vs_REAN.py -r $REANDIR -i $STAT_PROFILES_DIR -o Fig4.21/


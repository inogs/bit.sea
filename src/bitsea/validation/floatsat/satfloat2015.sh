# from qualification_quid.sh in ../deliverables

RUN=DA_FLOAT_SAT/Winter/RUN_SAT_01
export MASKFILE=/pico/scratch/userexternal/ateruzzi/$RUN/wrkdir/MODEL/meshmask.nc

#         INPUTDIR=/pico/scratch/userexternal/ateruzzi/DA_FLOAT_SAT/$RUN/wrkdir/MODEL/AVE_FREQ_1/
#STAT_PROFILES_DIR=/pico/scratch/userexternal/ateruzzi/eas_v12/eas_v19_3/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES

#mkdir -p Fig4.1/ Fig4.7 Fig4.9 Fig4.10 Fig4.17 Fig4.18
#COMMONS_PARAMS="-m $MASKFILE -o LayerMaps/  -l Plotlist_bio.xml -s 20150101 -e 20170101"

mkdir -p Fig4.2 Fig4.3/offshore Fig4.3/coast table4.1 table4.2
SAT_WEEKLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MULTISENSOR/DT/WEEKLY/
INPUT_AGGR_DIR=/pico/scratch/userexternal/ateruzzi/$RUN/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/
STAT_PROFILES_DIR1=/pico/scratch/userexternal/ateruzzi/$RUN/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/
OUTDIR=TMP_week/
python ScMYvalidation_plan.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c open_sea   -o export_data_ScMYValidation_plan_open_sea.pkl
python ScMYvalidation_plan.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c coast      -o export_data_ScMYValidation_plan_coast.pkl

exit 0

python plot_timeseries_2014_2016.py -i $STAT_PROFILES_DIR1 -m $MASKFILE -o export_data_ScMYValidation_plan_open_sea.pkl -c export_data_ScMYValidation_plan_coast.pkl -O ./Fig4.2/

# CHL-SURF-W-CLASS4-SIMG-BIAS-BASIN
# CHL-SURF-W-CLASS4-SIMG-RMSD-BASIN
python plot_timeseries_RMS.py -i export_data_ScMYValidation_plan_open_sea.pkl -o Fig4.3/offshore  # table4.1
python plot_timeseries_RMS.py -i export_data_ScMYValidation_plan_coast.pkl    -o Fig4.3/coast     # table4.2
cp Fig4.3/offshore/table4.1.dat table4.1
cp Fig4.3/coast/table4.1.dat    table4.2/table4.2.dat


mkdir -p Fig4.8
INTEGRALS_PPN=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v19_3/wrkdir/POSTPROC/output/AVE_FREQ_2/PPN/INTEGRALS/PPN/ # girato un aveScan ridotto solo per loro
python read_ppn_from_avescan_do_plot.py -c open_sea   -i $INTEGRALS_PPN -o Fig4.8


# Figure carbonatiche 1x1
# ALK-LAYER-Y-CLASS4-CLIM-RMSD
# ALK-LAYER-Y-CLASS4-CLIM-BIAS
# DIC-LAYER-Y-CLASS4-CLIM-RMS
# DIC-LAYER-Y-CLASS4-CLIM-BIAS

# ---------- executed elsewhere
# python ricostruzione_Integrals.py -i /pico/scratch/userexternal/gbolzon0/eas_v12/eas_v19_3/wrkdir/POSTPROC/output/AVE_FREQ_2/1x1/INTEGRALS/ -o 1x1/
# ----------------------------

# PH-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN in Fig4.20
mkdir -p Fig4.20/Feb Fig4.20/May Fig4.20/Aug Fig4.20/Nov Fig4.21/Feb Fig4.21/May  Fig4.21/Aug  Fig4.21/Nov
COMMONS_PARAMS="-i $INPUTDIR -m $MASKFILE -l Plotlist_bio.xml  -t mean "
python averager_and_plot_map.py $COMMONS_PARAMS  -v PH  -s 20160201 -e 20160301 -o Fig4.20/Feb
python averager_and_plot_map.py $COMMONS_PARAMS  -v PH  -s 20160501 -e 20160601 -o Fig4.20/May
python averager_and_plot_map.py $COMMONS_PARAMS  -v PH  -s 20160801 -e 20160901 -o Fig4.20/Aug
python averager_and_plot_map.py $COMMONS_PARAMS  -v PH  -s 20161101 -e 20161201 -o Fig4.20/Nov

# PCO-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN in Fig4.21
python averager_and_plot_map.py $COMMONS_PARAMS  -v pCO2  -s 20160201 -e 20160301 -o Fig4.21/Feb
python averager_and_plot_map.py $COMMONS_PARAMS  -v pCO2  -s 20160501 -e 20160601 -o Fig4.21/May
python averager_and_plot_map.py $COMMONS_PARAMS  -v pCO2  -s 20160801 -e 20160901 -o Fig4.21/Aug
python averager_and_plot_map.py $COMMONS_PARAMS  -v pCO2  -s 20161101 -e 20161201 -o Fig4.21/Nov


# BIOFLOATS SECTION: Hovmoeller plots, wmo trajectories and statistics per basin
# It is included the 2014 also.
# Figures 4.4a
mkdir -p Fig4.4a Fig4.4b Fig4.6 Fig4.12 Fig4.13 Fig4.15 tmp_nc table4.4 table4.8 
OUTDIR=Fig4.4a
#python Hov_flots+model.py -m $MASKFILE -o $OUTDIR
python Hov_flots+model_vars.py -m $MASKFILE -o $OUTDIR
mv $OUTDIR/*N3n*.png Fig4.12
mv $OUTDIR/*O2o*.png Fig4.15

# Figures 4.4b
NCDIR=tmp_nc
OUTDIR=Fig4.4b
python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -o $NCDIR
python SingleFloat_vs_Model_Stat_Timeseries_plotter.py -i $NCDIR -o $OUTDIR
mv $OUTDIR/N3n*.png Fig4.12
mv $OUTDIR/O2o*.png Fig4.15

# Figures 4.6 and 4.13 + tables 4.4 and 4.8
# CHL-PROF-D-CLASS4-PROF-CORR-BASIN
# NIT-PROF-D-CLASS4-PROF-CORR-BASIN
#  DO-PROF-D-CLASS4-PROF-CORR-BASIN 

OUTDIR=Fig4.6
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

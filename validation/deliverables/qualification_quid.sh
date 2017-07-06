# Descriptor of CMEMS-Med-biogeochemistry-ScQP-v1.3.docx

# QUID REANALYSIS
# SECTION 4?:
export MASKFILE=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/meshmask.nc

         INPUTDIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v19_2/wrkdir/MODEL/AVE_FREQ_2
   INPUT_AGGR_DIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v19_2/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP
STAT_PROFILES_DIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v19_2/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES
SAT_MONTHLY_DIR=

OUTDIR=spaghettiplots
mkdir -p $OUTDIR
python plot_layer_timeseries_on_profiles.py -i $STAT_PROFILES_DIR -m $MASKFILE -o $OUTDIR

mkdir -p Fig4.1/


COMMONS_PARAMS="-m $MASKFILE -o LayerMaps/  -l Plotlist_bio.xml -s 20150101 -e 20170101"

mkdir -p LayerMaps
python averager_and_plot_map.py -i $INPUT_AGGR_DIR  -v P_l  -t mean $COMMONS_PARAMS      # Fig4.1 CHL-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v N3n  -t mean $COMMONS_PARAMS      # FIg4.10 NIT-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v N1p  -t mean $COMMONS_PARAMS      # Fig4.9  PHOS-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v ppn  -t integral $COMMONS_PARAMS  # Fig4.7 per lo 0-200m

python averager_and_plot_map.py -i $INPUTDIR        -v ppn  -t mean $COMMONS_PARAMS  # NPP-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v O2o  -t mean $COMMONS_PARAMS  # DO-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUT_AGGR_DIR  -v P_c  -t mean $COMMONS_PARAMS  # PHYC-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v pH   -t mean $COMMONS_PARAMS  # PH-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v pCO2 -t mean $COMMONS_PARAMS  # PCO-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN


python sat_ave_and_plot.py      -i $SAT_MONTHLY_DIR -m $MASKFILE  -o Fig4.1/

mkdir -p fig4.2 fig4.3/offshore fig4.3/coast table4.1 table4.2
SAT_WEEKLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/WEEKLY_24_Friday
INPUT_AGGR_DIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_11/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP
python ScMYvalidation_plan.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c open_sea   -o export_data_ScMYValidation_plan_open_sea.pkl
python ScMYvalidation_plan.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c coast      -o export_data_ScMYValidation_plan_coast.pkl
python plot_timeseries.py -o export_data_ScMYValidation_plan_open_sea.pkl -c export_data_ScMYValidation_plan_coast.pkl -O ./fig4.2/

python plot_timeseries_RMS.py -i export_data_ScMYValidation_plan_open_sea.pkl -o fig4.3/offshore  # table4.1
python plot_timeseries_RMS.py -i export_data_ScMYValidation_plan_coast.pkl    -o fig4.3/coast     # table4.2
cp fig4.3/offshore/table4.1 table4.1
cp fig4.3/coast/table4.1    table4.2


mkdir -p Fig4.8
INTEGRALS_PPN=/gpfs/scratch/userexternal/gbolzon0/RA_COAST_ATM/wrkdir/POSTPROC/output/AVE_FREQ_2/ONLY_PPN/INTEGRALS/PPN/ # girato un aveScan ridotto solo per loro
python read_ppn_from_avescan_do_plot.py -c open_sea   -i $INTEGRALS_PPN -o Fig4.8


# Figure carbonatiche 1x1
# ALK-LAYER-Y-CLASS4-CLIM-RMSD
# ALK-LAYER-Y-CLASS4-CLIM-BIAS
# DIC-LAYER-Y-CLASS4-CLIM-RMS
# DIC-LAYER-Y-CLASS4-CLIM-BIAS
python ricostruzione_Integrals.py -i /gpfs/scratch/userexternal/gbolzon0/RA_COAST_02/wrkdir/POSTPROC/output/AVE_FREQ_2/1x1/INTEGRALS/ -o 1x1/

mkdir -p Fig4.19/Feb Fig1.19/May Fig4.19/Aug Fig4.19/Nov Fig4.20/Feb Fig4.20/May  Fig4.20/Aug  Fig4.20/Nov
python averager_and_plot_map.py -i $INPUTDIR  -v pH  -t mean -s 20160201 -e 20160301  -m $MASKFILE -o Fig4.19/Feb
python averager_and_plot_map.py -i $INPUTDIR  -v pH  -t mean -s 20160501 -e 20160601  -m $MASKFILE -o Fig4.19/May
python averager_and_plot_map.py -i $INPUTDIR  -v pH  -t mean -s 20160801 -e 20160901  -m $MASKFILE -o Fig4.19/Aug
python averager_and_plot_map.py -i $INPUTDIR  -v pH  -t mean -s 20161101 -e 20161201  -m $MASKFILE -o Fig4.19/Nov

python averager_and_plot_map.py -i $INPUTDIR  -v pCO2  -t mean -s 20160201 -e 20160301  -m $MASKFILE -o Fig4.20/Feb
python averager_and_plot_map.py -i $INPUTDIR  -v pCO2  -t mean -s 20160501 -e 20160601  -m $MASKFILE -o Fig4.20/May
python averager_and_plot_map.py -i $INPUTDIR  -v pCO2  -t mean -s 20160801 -e 20160901  -m $MASKFILE -o Fig4.20/Aug
python averager_and_plot_map.py -i $INPUTDIR  -v pCO2  -t mean -s 20161101 -e 20161201  -m $MASKFILE -o Fig4.20/Nov



# BIOFLOATS SECTION: Hovmoeller plots, wmo trajectories and statistics per basin
# Figures 4.4a
mkdir -p Fig4.4a Fig4.4b Fig4.5 Fig4.13 tmp_nc table4.4 table4.7
OUTDIR=Fig4.4a
python Hov_flots+model.py -m $MASKFILE -o $OUTDIR

# Figures 4.4b
NCDIR=tmp_nc
OUTDIR=Fig4.4b
python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -o $NCDIR
python SingleFloat_vs_Model_Stat_Timeseries_plotter.py -i $NCDIR -o $OUTDIR
mv $OUTDIR/N3n*.png Fig4.12

# Figures 4.5 and 4.13 + tables 4.4 and 4.7
# CHL-PROF-D-CLASS4-PROF-CORR-BASIN
# NIT-PROF-D-CLASS4-PROF-CORR-BASIN
#  DO-PROF-D-CLASS4-PROF-CORR-BASIN 

OUTDIR=Fig4.5
python BASIN_Float_vs_Model_Stat_Timeseries_monthly.py -m $MASKFILE -o $NCDIR
python BASIN_Float_vs_Model_Stat_Timeseries_monthly_plotter.py -m $MASKFILE -i $NCDIR -o $OUTDIR
mv $OUTDIR/N3n*.png Fig4.13

cp $OUTDIR/P_l_tab_statistics_SHORT.txt table4.4/ 
cp $OUTDIR/N3n_tab_statistics_SHORT.txt table4.7/


# BIOFLOATS SECTION: statistics on layers
# CHL-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN
# NIT-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN
#  DO-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN


OUTFIGDIR=Floats_bias_rmse_Timeseries     # 8layer x 7sub x 3var = 168 png files
TABLE_DIR=Floats_bias_rmse_tables         #: 2stats x 3var        = 6 txt files, TABLE.O2o_BIAS.txt  with time average for each layer,sub
mkdir -p $OUTFIGDIR $TABLE_DIR table4.3/ table4.8/ table4.11/
python biofloats_ms.py  -m $MASKFILE -o float_bias_rmse.nc
python biofloats_ms_plotter.py -i float_bias_rmse.nc -f $OUTFIGDIR -t $TABLE_DIR
cp $TABLE_DIR/P_l_BIAS.txt $TABLE_DIR/P_l_RMSE.txt table4.3/
cp $TABLE_DIR/N3n_BIAS.txt $TABLE_DIR/N3n_RMSE.txt table4.8/
cp $TABLE_DIR/O2o_BIAS.txt $TABLE_DIR/O2o_RMSE.txt table4.11/
cp $OUTFIGDIR/*P_l* fig4.5
cp $OUTFIGDIR/*N3n* fig4.14
cp $OUTFIGDIR/*O2o* fig4.15


#########################   static dataset climatology section ###################################
# Figures 4.11 and 4.18
mkdir -p sim_vs_clim_profiles/ Fig4.11 Fig4.18
python simulation_vs_clim.py -i $STAT_PROFILES_DIR -o sim_vs_clim_profiles/ -s 20150101 -e 20170101 -m $MASKFILE
cp sim_vs_clim_profiles/Fig_4.11*png Fig.4.11
cp sim_vs_clim_profiles/Fig_4.18*png Fig.4.18

DIR=static_clim
mkdir -p $DIR table4.5 table4.6 table4.9/ table4.10 table4.12/ table4.13
# -------------------------------------------------------------------------

python static_clim_validation.py -i $STAT_PROFILES_DIR -o $DIR -m $MASKFILE -s 20150101 -e 20170101

# -------------------------------------------------------------------------
cp $DIR/N1p-LAYER-Y-CLASS4-CLIM.txt $DIR/N3n-LAYER-Y-CLASS4-CLIM.txt table4.5/
cp $DIR/O2o-LAYER-Y-CLASS4-CLIM.txt                                  table4.9/
cp $DIR/Ac-LAYER-Y-CLASS4-CLIM.txt $DIR/DIC-LAYER-Y-CLASS4-CLIM.txt  table4.12/

cp $DIR/N1p-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt $DIR/N3n-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt table4.6
cp $DIR/N1p-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt                                            table4.10
cp $DIR/Ac-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt  $DIR/DIC-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt table4.13

# PHO-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR
# NIT-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR
#  DO-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR
# ALK-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR
# DIC-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR  calculated on 8 layers

# ALK-PROF-Y-CLASS4-CLIM-CORR-BASIN
# DIC-PROF-Y-CLASS4-CLIM-CORR-BASIN      calculated on 14 layers


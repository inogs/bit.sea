# Descriptor of CMEMS-Med-biogeochemistry-ScQP-v1.3.docx

# QUID REANALYSIS
# SECTION 4?:
export MASKFILE=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/meshmask.nc
SAT_MONTHLY_DIR=

INPUTDIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/AVE_FREQ_2
INPUT_AGGR_DIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP



COMMONS_PARAMS="-m $MASKFILE -o LayerMaps/  -l Plotlist_bio.xml -s 20140101 -e 20180101"

mkdir -p LayerMaps
python averager_and_plot_map.py -i $INPUT_AGGR_DIR  -v P_l  -t mean $COMMONS_PARAMS  # CHL-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v N3n  -t mean $COMMONS_PARAMS  # NIT-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v N1p  -t mean $COMMONS_PARAMS  # PHOS-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v ppn  -t mean $COMMONS_PARAMS  # NPP-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v O2o  -t mean $COMMONS_PARAMS  # DO-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUT_AGGR_DIR  -v P_c  -t mean $COMMONS_PARAMS  # PHYC-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v pH   -t mean $COMMONS_PARAMS  # PH-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v pCO2 -t mean $COMMONS_PARAMS  # PCO-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v ppn  -t integral $COMMONS_PARAMS  # per lo 0-200m

SAT_WEEKLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MULTISENSOR_1km/WEEKLY_24_Friday
INPUT_AGGR_DIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_11/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP
python ScMYvalidation_plan.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c open_sea   -o export_data_ScMYValidation_plan_open_sea.pkl
python ScMYvalidation_plan.py -s $SAT_WEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c coast      -o export_data_ScMYValidation_plan_coast.pkl
python plot_timeseries.py -o export_data_ScMYValidation_plan_open_sea.pkl -c export_data_ScMYValidation_plan_coast.pkl -O ./fig4.2/


STAT_PROFILES_DIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_11/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES
OUTDIR=spaghettiplots
mkdir -p $OUTDIR
python plot_layer_timeseries_on_profiles.py -i $STAT_PROFILES_DIR -m $MASKFILE -o $OUTDIR

#bioflots section
OUTDIR=Timeseries_rmse_floats
mkdir -p $OUTDIR
python biofloats_ms.py  -m $MASKFILE -o float_bias_rmse.nc  #  CHL-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN
python biofloats_ms_plotter.py -i float_bias_rmse.nc -o $OUTDIR


NCDIR=/pico/scratch/userexternal/lfeudale/validation/work/output/
OUTDIR=/pico/scratch/userexternal/lfeudale/validation/work/output/PNG/
python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -o $NCDIR
python SingleFloat_vs_Model_Stat_Timeseries_plotter.py -i $NCDIR -o $OUTDIR





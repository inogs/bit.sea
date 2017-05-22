# Descriptor of CMEMS-Med-biogeochemistry-ScQP-v1.3.docx

# QUID REANALYSIS
# SECTION 4?:
export MASKFILE=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/meshmask.nc
SAT_MONTHLY_DIR=

INPUTDIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/AVE_FREQ_2
INPUT_AGGR_DIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP



COMMONS_PARAMS="-m $MASKFILE -o LayerMaps/ -t mean -l Plotlist_bio.xml -s 20140101 -e 20180101"

mkdir -p LayerMaps
python averager_and_plot_map.py -i $INPUT_AGGR_DIR  -v P_l  $COMMONS_PARAMS  # CHL-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v N3n  $COMMONS_PARAMS  # NIT-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN 
python averager_and_plot_map.py -i $INPUTDIR        -v N1p  $COMMONS_PARAMS  # PHOS-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v ppn  $COMMONS_PARAMS  # NPP-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v O2o  $COMMONS_PARAMS  # DO-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUT_AGGR_DIR  -v P_c  $COMMONS_PARAMS  # PHYC-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v pH   $COMMONS_PARAMS  # PH-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
python averager_and_plot_map.py -i $INPUTDIR        -v pCO2 $COMMONS_PARAMS  # PCO-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
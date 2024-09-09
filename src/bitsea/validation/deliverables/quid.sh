# Descriptor of CMEMS-Med-QUID-006-008-V2-V1.0.docx

# QUID REANALYSIS
# SECTION 4?:
#export MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc
export MASKFILE=/gpfs/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc
#SAT_MONTHLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI/MONTHLY_V4/
SAT_MONTHLY_DIR=/gpfs/scratch/userexternal/ateruzzi/SAT_forRA_COAST/TIMESER_SATMONTHLY_RACOAST/

#INPUTDIR=/pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/
#INPUTDIR=/gpfs/scratch/userexternal/gbolzon0/RA_COAST/wrkdir/MODEL/AVE_FREQ_2/
#INPUT_AGGR_DIR=/gpfs/scratch/userexternal/gbolzon0/RA_COAST/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP
INPUT_AGGR_DIR=/gpfs/scratch/userexternal/ateruzzi/RA_COAST_P_l/

mkdir -p fig4.2  table4.3  Fig4.1  Fig4.6 Fig4.5 Fig4.4 Fig4.7 Fig4.8 Fig4.9 Fig4.10 Fig4.11 Fig4.12
mkdir -p table4.3 table4.4 table4.5
python ScMYvalidation_plan.py -s $SAT_MONTHLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c open_sea   -o export_data_ScMYValidation_plan_open_sea.pkl
python ScMYvalidation_plan.py -s $SAT_MONTHLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c coast      -o export_data_ScMYValidation_plan_coast.pkl
python ScMYvalidation_plan.py -s $SAT_MONTHLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE -c everywhere -o export_data_ScMYValidation_plan_everywhere.pkl

python plot_timeseries.py -o export_data_ScMYValidation_plan_open_sea.pkl -c export_data_ScMYValidation_plan_coast.pkl -O ./fig4.2/

# figure 4.3 and table 4.1
mkdir -p fig4.3/offshore fig4.3/coast
python plot_timeseries_RMS.py -i export_data_ScMYValidation_plan_open_sea.pkl -o fig4.3/offshore
python plot_timeseries_RMS.py -i export_data_ScMYValidation_plan_coast.pkl    -o fig4.3/coast

python MYvalidation_statics.py -m $MASKFILE -o export_data_ScMYValidation_plan_statics.pkl
# table 4.3 and 4.4
python reader_statics.py -o ./table4.3 # phosphate nitrate o2


python averager_and_plot_map.py -i $INPUT_AGGR_DIR  -m $MASKFILE  -o Fig4.1/ -v P_l -t mean -l Plotlist_bio.xml -s 19990101 -e 20171231
python sat_ave_and_plot.py      -i $SAT_MONTHLY_DIR -m $MASKFILE  -o Fig4.1/

python averager_and_plot_map.py -i $INPUTDIR -m $MASKFILE -o Fig4.6/ -v N1p -t mean  -l nut_layerlist
python averager_and_plot_map.py -i $INPUTDIR -m $MASKFILE -o Fig4.6/ -v N3n -t mean  -l nut_layerlist
#questi poi sono da modificare in maniera da avere questi limiti :
#N1p sempre [0,0.15], ma dai 100-150 in poi diventa [0,0.35] ,
#N3n sempre [0,4]   , da dai 60-100  in poi diventa [0, 10]


mkdir FigureP_c
python averager_and_plot_map.py -i $INPUT_AGGR_DIR -m $MASKFILE -o FigureP_c -v P_c -t mean  -l nut_layerlist

# Figure IV.5
INTEGRALS_PPN=/gpfs/scratch/userexternal/gbolzon0/RA_COAST_ATM/wrkdir/POSTPROC/output/AVE_FREQ_2/11_sub/INTEGRALS/PPN/
INTEGRALS_PPN=/gpfs/scratch/userexternal/gbolzon0/RA_COAST_ATM/wrkdir/POSTPROC/output/AVE_FREQ_2/ONLY_PPN/INTEGRALS/PPN/ # girato un aveScan ridotto solo per loro

mkdir -p Fig4.5/coast Fig4.5/offshore Fig4.5/everywhere
python read_ppn_from_avescan_do_plot.py -c coast      -i $INTEGRALS_PPN -o Fig4.5/coast
python read_ppn_from_avescan_do_plot.py -c open_sea   -i $INTEGRALS_PPN -o Fig4.5/offshore
python read_ppn_from_avescan_do_plot.py -c everywhere -i $INTEGRALS_PPN -o Fig4.5/everywhere

python averager_and_plot_map.py -i $INPUTDIR -m $MASKFILE  -o Fig4.4/ -v ppn -t integral -l ppn_layerlist


#Figure IV.7 - density PHOSPHATE
python density_plots.py     -M $MASKFILE -o Fig4.7/  -v N1p -m 0
python vertical_profiles.py -m $MASKFILE -o Fig4.8/  -v N1p | grep corr > ./table4.3/table.4.3_corr.dat

python density_plots.py     -M $MASKFILE -o Fig4.9/  -v N3n -m 0
python vertical_profiles.py -m $MASKFILE -o Fig4.10/ -v N3n | grep corr > ./table4.4/table4.4_corr.dat

python density_plots.py     -M $MASKFILE -o Fig4.11/ -v O2o
python vertical_profiles.py -m $MASKFILE -o Fig4.12/ -v O2o | grep corr > ./table4.5/table4.5_corr.dat

# Figure carbonatiche
# Fig 4.13, 4.14, per AC_ e DIC
# readMAP1x1_13layer_do_CFR_carbsys.m usa questi
# /pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/POSTPROC/output/AVE_FREQ_2/1x1/MAPS/MAP1x1_13lev_' + varname +'.nc
# Generazione (in una bin_13lev_1x1)
# aveScan.py con maskload preso da carbonatics/maskload_1x1.py
#  in getAllStatistics:
#   - solo integrali
#   - dopo m = SUB(sub) uscire così [if m.sum() ==0 : continue ] perché ci sono troppi sottobacini di terra
#  in Vol_Integrals
#   - chiamare CoreStatistics_noSort:
python ricostruzione_Integrals.py -i /gpfs/scratch/userexternal/gbolzon0/RA_COAST_02/wrkdir/POSTPROC/output/AVE_FREQ_2/1x1/INTEGRALS/ -o 1x1/



# Fig 4.15, 4.16
# readQUADRATI4x4_PROFILI_do_plotPROFILI_monovariate.m
# che usa PROF_18aree_${VAR}.nc

# Generazione (in una bin_4x4)
# aveScan.py con maskload preso da carbonatics/maskload_4x4.py
#  in getAllStatistics:
#   - solo profili

for var in Ac DIC O3c O3h pCO2 PH P_l; do
    python ricostruzione_profili.py -i /gpfs/scratch/userexternal/gbolzon0/RA_COAST_02/wrkdir/POSTPROC/output/AVE_FREQ_2/4x4/STAT_PROFILES -m $MASKFILE -v $var
done


mkdir ./Fig4.17
python seasonal_plot_map.py -i $INPUTDIR -o Fig4.17/ -m 2 -v PH -t mean -M $MASKFILE
python seasonal_plot_map.py -i $INPUTDIR -o Fig4.17/ -m 5 -v PH -t mean -M $MASKFILE
python seasonal_plot_map.py -i $INPUTDIR -o Fig4.17/ -m 8 -v PH -t mean -M $MASKFILE
python seasonal_plot_map.py -i $INPUTDIR -o Fig4.17/ -m 11 -v PH -t mean -M $MASKFILE

python seasonal_plot_map.py -i $INPUTDIR -o Fig4.17/ -m 2 -v pCO2 -t mean -M $MASKFILE
python seasonal_plot_map.py -i $INPUTDIR -o Fig4.17/ -m 5 -v pCO2 -t mean -M $MASKFILE
python seasonal_plot_map.py -i $INPUTDIR -o Fig4.17/ -m 8 -v pCO2 -t mean -M $MASKFILE
python seasonal_plot_map.py -i $INPUTDIR -o Fig4.17/ -m 11 -v pCO2 -t mean -M $MASKFILE


#--------------------------------------------------
#QUID ANALYSIS AND FORECAST
Figure IV.5.
Figure IV.6
/pico/home/userexternal/gcossari/bit.sea/calibration_bioFloat.py

Figure IV.9.
Figure IV.10.
/pico/home/userexternal/gcossari/bit.sea/calibration_mooring.py

Figure IV.11.
Figure IV.12.
Su pico
/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/readVAR_doPROFILI_18aree.py
/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/readVAR_doMAP1DEG_13LAYERs.py
e
read_mask1x1_npy_do_netcdf.py


poi i netcdf sono usati in locale questi due matlab:
readMAP1x1_13layer_do_CO2sys.m
readQUADRATI4x4_PROFILI_do_plotPROFILI.m

# Validatio of coastal assimilation runs

export MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc
SAT_WEEKLY_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI/NEW_20161702/WEEKLY_V4/

RUN=DA_COAST_01
MOD_WEEKLY_DIR=/pico/scratch/userexternal/ateruzzi/$RUN/wrkdir/POSTPROC/output/DA__FREQ_1/TMP/

OUTDIR=/pico/scratch/userexternal/ateruzzi/ELAB_DA_COAST/OUTPUTvalidation/$RUN
mkdir $OUTDIR

# Time series
mkdir $OUTDIR/TIMESER/
python ScMYvalidation_plan.py -o export_data_ScMYValidation_plan.pkl -s $SAT_WEEKLY_DIR -m $MOD_WEEKLY_DIR

python plot_timeseries.py -i export_data_ScMYValidation_plan.pkl -o $OUTDIR/TIMESER/


# RMS and BOAS time series
mkdir $OUTDIR/BIAS_RMS/
python plot_timeseries_RMS.py -i export_data_ScMYValidation_plan.pkl -o $OUTDIR/BIAS_RMS/


mkdir $OUTDIR/table4.3/
python MYvalidation_statics.py -o export_data_ScMYValidation_plan_statics.pkl
# table 4.3 and 4.4
python reader_statics.py -o $OUTDIR/table4.3 # phosphate nitrate o2


#Figure 4.1
mkdir $OUTDIR/Fig4.1

python averager_and_plot_map.py -i $MOD_MONTHLY_DIR       -o $OUTDIR/Fig4.1/ -v P_i -t mean # Time interval hard coded!
# e per avere solo il primo livello?
# aggiungo argomento mese
python sat_ave_and_plot.py      -i $SAT_MONTHLY_DIR       -o $OUTDIR/Fig4.1/

mkdir $OUTDIR/Fig4.6
python averager_and_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.6/ -v N1p -t mean  --top 0 --bottom  50 --mapdepthfilter  50.0
python averager_and_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.6/ -v N1p -t mean  --top 0 --bottom 150 --mapdepthfilter 150.0

python averager_and_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.6/ -v N3n -t mean  --top 0 --bottom   50 --mapdepthfilter  50.0
python averager_and_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.6/ -v N3n -t mean  --top 0 --bottom  150 --mapdepthfilter 150.0

# Figure IV.5
python read_ppn_from_avescan_do_plot.py -i /pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/POSTPROC/output/AVE_FREQ_2/only_ppn/INTEGRALS/PPN/ -o $OUTDIR/Fig4.5/

mkdir $OUTDIR/Fig4.4
python averager_and_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.4/ -v ppn -t integral --top 0 --bottom 200 --mapdepthfilter 150.0

mkdir $OUTDIR/table4.4 $OUTDIR/table4.5
#Figure IV.7 - density PHOSPHATE
python density_plots.py     -o $OUTDIR/Fig4.7/  -v N1p -m 0
python vertical_profiles.py -o $OUTDIR/Fig4.8/  -v N1p > ./table4.3/table.4.3_corr.dat

python density_plots.py     -o $OUTDIR/Fig4.9/  -v N3n -m 0
python vertical_profiles.py -o $OUTDIR/Fig4.10/ -v N3n > ./table4.4/table4.4_corr.dat

python density_plots.py     -o $OUTDIR/Fig4.11/ -v O2o
python vertical_profiles.py -o $OUTDIR/Fig4.12/ -v O2o > ./table4.5/table4.5_corr.dat

# Figure carbonatiche
# Fig 4.13, 4.14, per AC_ e DIC
# readMAP1x1_13layer_do_CFR_carbsys.m usa questi
# /pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/POSTPROC/output/AVE_FREQ_2/1x1/MAPS/MAP1x1_13lev_' + varname +'.nc
# che vengono generati da ricostruzione_Integrals.py di opa_rea/chain/postproc
# che usa un maskload a parte e un aveScan che fa solo integrali, definiti qui
# /pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/POSTPROC/bin_13lev_1x1

# Fig 4.15, 4.16
# readQUADRATI4x4_PROFILI_do_plotPROFILI_monovariate.m
# che legge da qui
# /pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/POSTPROC/bin_4x4
# i files PROF_18aree_${VAR}.nc, generati da ricostruzione_profili.py

mkdir $OUTDIR/Fig4.17
python seasonal_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.17/ -m 2 -v PH -t mean
python seasonal_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.17/ -m 5 -v PH -t mean
python seasonal_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.17/ -m 8 -v PH -t mean
python seasonal_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.17/ -m 11 -v PH -t mean

python seasonal_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.17/ -m 2 -v pCO2 -t mean
python seasonal_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.17/ -m 5 -v pCO2 -t mean
python seasonal_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.17/ -m 8 -v pCO2 -t mean
python seasonal_plot_map.py -i $MOD_MONTHLY_DIR -o $OUTDIR/Fig4.17/ -m 11 -v pCO2 -t mean


--------------------------------------------------
QUID ANALYSIS AND FORECAST
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

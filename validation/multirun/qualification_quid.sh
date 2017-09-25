#! /bin/bash



mkdir -p Fig4.3/offshore Fig4.3/coast
V2_PKL=/pico/scratch/userexternal/lfeudale/validation/V2C/bit.sea/validation/deliverables/export_data_ScMYValidation_plan_open_sea.pkl
V3_PKL=/pico/scratch/userexternal/lfeudale/validation/eas_v12/eas_v20_1/bit.sea/validation/deliverables/export_data_ScMYValidation_plan_open_sea.pkl
python plot_timeseries_RMS_2.py -i1 $V3_PKL -i2 $V2_PKL    -o Fig4.3/offshore  # table4.1




V2_PKL=/pico/scratch/userexternal/lfeudale/validation/V2C/bit.sea/validation/deliverables/export_data_ScMYValidation_plan_coast.pkl
V3_PKL=/pico/scratch/userexternal/lfeudale/validation/eas_v12/eas_v20_1/bit.sea/validation/deliverables/export_data_ScMYValidation_plan_coast.pkl
python plot_timeseries_RMS_2.py -i1 $V3_PKL -i2 $V2_PKL   -o Fig4.3/coast     # table4.2


OUTFIGDIR=Floats_bias_rmse_Timeseries     # 8layer x 7sub x 3var = 168 png files
mkdir -p $OUTFIGDIR Fig4.5/ Fig4.14/ Fig4.16/

V2_NC=/pico/scratch/userexternal/lfeudale/validation/V2C/bit.sea/validation/deliverables/float_bias_rmse.nc
#V3_NC=/pico/scratch/userexternal/lfeudale/validation/eas_v12/eas_v20_1/bit.sea/validation/deliverables/float_bias_rmse.nc
V3_NC=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v20_1/wrkdir/POSTPROC/bit.sea/validation/deliverables/float_bias_rmse.nc

python biofloats_ms_plotter.py -i1 $V3_NC -i2 $V2_NC -f $OUTFIGDIR
cp $OUTFIGDIR/*P_l* Fig4.5
cp $OUTFIGDIR/*N3n* Fig4.14
cp $OUTFIGDIR/*O2o* Fig4.16



export MASKFILE=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v19_3/wrkdir/MODEL/meshmask.nc
STAT_PROFILES_MONTHLY_DIR_V3=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v20_1/wrkdir/POSTPROC/output/MONTHLY/STAT_PROFILES/
STAT_PROFILES_MONTHLY_DIR_V2=/pico/scratch/userexternal/gbolzon0/V2C/MONTHLY/STAT_PROFILES/
MASK_V2=/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc
mkdir -p sim_vs_clim_profiles/ Fig4.11 Fig4.19
python simulation_vs_clim.py -i1 $STAT_PROFILES_MONTHLY_DIR_V3 -m1 $MASKFILE -i2 $STAT_PROFILES_MONTHLY_DIR_V2 -m2 $MASK_V2   -o sim_vs_clim_profiles/ -s 20150101 -e 20170101 

cp sim_vs_clim_profiles/Fig_4.11*png Fig4.11
cp sim_vs_clim_profiles/Fig_4.19*png Fig4.19
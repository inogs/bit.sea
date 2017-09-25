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

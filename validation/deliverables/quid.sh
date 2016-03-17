# Descriptor of CMEMS-Med-QUID-006-008-V2-V1.0.docx

# QUID REANALYSIS
# SECTION 4?:
export MASKFILE=/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc

mkdir ./fig4.2/
python ScMYvalidation_plan.py -o export_data_ScMYValidation_plan.pkl -s /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI/NEW_20161702/MONTHLY_V4/

# figure 4.2
mkdir ./fig4.2/
python plot_timeseries.py -i export_data_ScMYValidation_plan.pkl -o ./fig4.2/


# figure 4.3 and table 4.1
mkdir ./fig4.3/
python plot_timeseries_RMS.py -i export_data_ScMYValidation_plan.pkl -o ./fig4.3/


mkdir ./table4.3/
python MYvalidation_statics.py -o export_data_ScMYValidation_plan_statics.pkl
# table 4.3 and 4.4
python reader_statics.py -o ./table4.3 # phosphate nitrate o2




# Figure IV.5
python read_ppn_from_avescan_do_plot.py -i /pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/POSTPROC/output/AVE_FREQ_2/only_ppn/INTEGRALS/PPN/ -o Fig4.5/
#Figure IV.4.
#Figure IV.1
#Figure IV.6.
/pico/home/userexternal/gcossari/bit.sea/averager_and_plot_map.py
MA CON DEFINIZIONI INTERNE DIVERSE


#Figure IV.7 - density PHOSPHATE
python density_plots.py     -o Fig4.7/  -v N1p
python vertical_profiles.py -o Fig4.8/  -v N1p

python density_plots.py     -o Fig4.9/  -v N3n
python vertical_profiles.py -o Fig4.10/ -v N3n

python density_plots.py     -o Fig4.11/ -v O2o
python vertical_profiles.py -o Fig4.12/ -v O2o

Figure carbonatiche…
devo verificare dove sono

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

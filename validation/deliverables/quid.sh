# Descriptor of CMEMS-Med-QUID-006-008-V2-V1.0.docx

# QUID REANALYSIS
# SECTION 4?:
python ScMYvalidation_plan.py -o export_data_ScMYValidation_plan.pkl -s /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI/NEW_20161702/MONTHLY_V4/

# figure 4.2
plot_timeseries.py -i export_data_ScMYValidation_plan.pkl -o dir1/


# figure 4.3
 plot_timeseries_RMS.py -i export_data_ScMYValidation_plan.pkl -o dir2/

# table 4.1
   plot_timeseries_RMS.py



MYValidation_plan_statics.py # writes export_data_ScMYValidation_plan_statics.pkl
# table 4.3
reader_statics.py# phosphate
# table 4.4
reader_statics.py# Nitrate



# Figure IV.5
/pico/home/userexternal/gcossari/bit.sea/read_ppn_from_avescan_do_plot.py 

#Figure IV.4.
#Figure IV.1
#Figure IV.6. 
/pico/home/userexternal/gcossari/bit.sea/averager_and_plot_map.py 
MA CON DEFINIZIONI INTERNE DIVERSE

Figure IV.7
Figure IV.8.
Figure IV.9.
Figure IV.10.
Figure IV.11. 
Figure IV.12. 
fatte da stefano

Figure carboniche…
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



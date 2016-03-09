# Descriptor of CMEMS-Med-QUID-006-008-V2-V1.0.docx


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

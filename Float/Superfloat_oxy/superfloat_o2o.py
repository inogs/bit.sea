import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates ave files for aveScan profiler in chain validation
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--time_start','-Tst',
                                type = str,
                                required = True,
                                help = 'start time serie')

    parser.add_argument(   '--time_end','-Tend',
                                type = str,
                                required = True,
                                help = 'end time serie')

    parser.add_argument(   '--day','-Tday',
                                type = str,
                                required = True,
                                help = 'selected day')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = ''' output dir for .nc 
                                           ''')
    parser.add_argument(   '--variable', '-v',
                                type = str,
                                default = None,
                                required = True,
                                help = '''model variable''')

    return parser.parse_args()

args = argument()
import numpy as np
import pandas as pd 
from commons import timerequestors
from instruments import superfloat
from instruments import superfloat as bio_float
from instruments.var_conversions import FLOATVARS
import basins.OGS as OGS
from commons.utils import addsep
import warnings
warnings.filterwarnings('ignore')
import CORIOLIS_checks 
import sys
from commons_ import col_to_dt
import TREND_ANALYSIS
import superfloat_oxy_function
import glob
import os

DATE_start  = args.time_start
DATE_end    = args.time_end
DATE_DAY    = args.day
OUTDIR      = addsep(args.outdir)
varmod      = args.variable

#DATE_start='20140209'
#DATE_end='20170210'
#DATE_DAY='20170209'
#OUTDIR='/g100_scratch/userexternal/camadio0/SUPERFLOAT_old/'
#varmod='O2o'

# INPUT
THRES      = 20  # threshold to interpolate data at 600m using  580-620 m layer
VARNAME='DOXY'
LIST_DEPTH  = [600,800]
OUT_META ='OUTPUTS/'


# Get argo time series for 3 years
TI              = timerequestors.TimeInterval(starttime=DATE_start, endtime=DATE_end, dateformat='%Y%m%d')
Profilelist     = superfloat.FloatSelector(FLOATVARS[varmod],TI, OGS.med)
Profilelist     = superfloat_oxy_function.remove_bad_sensors(Profilelist , varmod)
# Get argo list for the selected day (DATE_DAY)
TI_1            = timerequestors.TimeInterval(starttime=DATE_DAY , endtime=DATE_end, dateformat='%Y%m%d')
Profilelist_day = superfloat.FloatSelector(FLOATVARS[varmod],TI_1, OGS.med)
Profilelist_day = superfloat_oxy_function.remove_bad_sensors(Profilelist_day , varmod)
WMO_LIST        = superfloat.get_wmo_list(Profilelist_day)

# if no data in day exit
print ('________________________ START _______________________________')
CORIOLIS_checks.check_data(WMO_LIST, DATE_DAY)

# define the dimension of file to fill 
dim = 0
for WMO in WMO_LIST:
    wmo_timeseries = superfloat.filter_by_wmo(Profilelist,WMO)
    dim            = np.sum([dim, len(wmo_timeseries)*len(LIST_DEPTH)])

# interpolated at 600m and 800m --> df filled 
COUNT=0
columnlist = ['ID','time','lat','lon','name','Type','VAR', 'Depth']
df = pd.DataFrame(index=np.arange(0,dim), columns=columnlist)
for DEPTH in LIST_DEPTH:
    for WMO in WMO_LIST:
        wmo_timeseries = superfloat.filter_by_wmo(Profilelist,WMO)
        for profile in wmo_timeseries:
            Pres, Profile, Qc = profile.read(FLOATVARS[varmod])
            IDX = np.where((Pres >= DEPTH-THRES  ) & ( Pres <= DEPTH+THRES))
            lst = [profile.ID(), profile.time, profile.lat,profile.lon, profile.name(),
                  (profile._my_float.origin(VARNAME)).status_var, np.nan, np.nan]
            df.loc[COUNT] = pd.Series(lst , columnlist)
            if np.array(IDX).size ==0:
               df.Depth[COUNT] = DEPTH
            else:
               df.Depth[COUNT] = DEPTH
               df.VAR[COUNT]   = np.interp(DEPTH , Pres[IDX], (Profile[IDX]))
            COUNT+=1

dim         = len(WMO_LIST) * len(LIST_DEPTH)
COLUMNS     = ['WMO','Depth','DURATION','min_date','max_date','Theil-Sen','RANSAC']
df_report   = pd.DataFrame(index=np.arange(0,dim), columns=COLUMNS)
from TREND_ANALYSIS import trend_conditions as TD
COUNT=0
TI_3            = timerequestors.TimeInterval(starttime=DATE_DAY, endtime=DATE_DAY, dateformat='%Y%m%d')

for DEPTH in LIST_DEPTH:
    for WMO in WMO_LIST:
        tmp = df[(df.Depth == DEPTH) & (df.name == WMO)]
        tmp.index = np.arange(0,len(tmp.index))
        days, min_d , max_d = CORIOLIS_checks.lenght_timeseries(tmp, 'time')
        if days > 0:
           Bool = CORIOLIS_checks.nans_check(tmp, 'VAR')
           tmp.dropna(inplace=True) # if no statistical estimators functin doesn work 
           lst = TD(WMO,days,Bool,DEPTH,min_d, max_d , TI_3, tmp)
        else: 
           # Argo 1st day 
           lst   = [WMO ,DEPTH, days , np.nan , np.nan , np.nan , np.nan ]
        df_report.iloc[COUNT,:] = pd.Series(lst, df_report.columns)
        COUNT+=1

df_report['TREND_TIME_SERIES'] ,df_report['TREND_per_YEAR'] ,df_report['DRIFT_CODE'] = np.nan , np.nan , np.nan 

from TREND_ANALYSIS import sign_analysis
for WMO in WMO_LIST:
    tmp  = df_report.loc[df_report.WMO==WMO]
    A    = np.append(np.array( tmp['Theil-Sen']), np.array( tmp['RANSAC']))
    if tmp.DURATION.any()==0:
        pass
    else:
        Bool = sign_analysis(A)
        df_report =  TREND_ANALYSIS.drift_coding(WMO, Bool, tmp, df_report)

Saving_report=False
if Saving_report:
   df_report.to_csv(OUT_META+DATE_start+'report_600_800.csv')

plot_time_series=False
import commons_
if plot_time_series:
   commons_.time_serie_plot(LIST_DEPTH, WMO_LIST, df , df_report) 
# EMODNET CHECK ONLY AT DEPTH =600
DEPTH = 600
df        =  df[df.Depth==DEPTH]
df_report =  df_report[df_report.Depth == DEPTH]
df        =  col_to_dt(df,'time')
dfday     =  df[df.date==TI_3.start_time]

# emodnet climatology iNIT
from basins.region import Region, Rectangle
from commons_ import cross_Med_basins
from commons_ import drop_rows_df

df_clim = pd.read_csv('EMODNET_climatology.csv',index_col=0)
df_cstd = pd.read_csv('EMODNET_stdev.csv',index_col=0)
COLINDEX=0

df_report['OFFSET'], df_report['basin'], df_report['EMODNET'], df_report['STdev*2'], df_report['Black_list'], df_report['filename'] = np.nan , np.nan , np.nan , np.nan, False , np.nan

for WMO in WMO_LIST:
    tmp_0    = dfday.loc[dfday.name == WMO]
    tmp_meta = df_report[ df_report.WMO==WMO]
    profile =  superfloat.filter_by_wmo(Profilelist_day ,str(WMO))
    if len(profile)>1:
        print (str(len(profile)) + ' profile available for ' +WMO +' in date '+ DATE_DAY )
        Saving_hr_csv=True
        if Saving_hr_csv:
           OUTPATH_NAME =  OUT_META + "High_time_freq_argo.csv" 
           indexlenght  = len(profile)
           columns_list=['WMO', 'DATE_DAY']
           values_list = [int(WMO), DATE_DAY ]
           commons_.save_report(OUTPATH_NAME,indexlenght,columns_list, values_list )
    for NPROF in profile:
        Pres, Profile, Qc = NPROF.read(FLOATVARS[varmod])
        ID_prof           = NPROF.ID()
        tmp               = tmp_0[tmp_0.ID==ID_prof]
        if tmp.VAR.isnull().values.any() or tmp_0.empty:
           print('no data for -->'  + WMO)
           outname  = superfloat_oxy_function.get_outfile(NPROF, OUTDIR)
           df_report.loc[ df_report.WMO == WMO , 'filename'] =  outname
           Condition_to_write=True
           if Condition_to_write:
              SAVEDIR = OUTDIR + WMO
              if not os.path.exists(SAVEDIR):
                 os.makedirs(SAVEDIR)
              import superfloat_generator
              writing_mode=superfloat_generator.writing_mode(outname)
              condition_to_write = superfloat_generator.exist_valid_variable('DOXY',outname)
              metadata = NPROF._my_float
              superfloat_oxy_function.dump_oxygen_file(outname, NPROF , Pres, Profile, Qc, metadata, mode=writing_mode)
        else:
           ARGO     = Rectangle(np.float(tmp.lon) , np.float(tmp.lon) , np.float(tmp.lat) , np.float(tmp.lat))
           NAME_BASIN, BORDER_BASIN  = cross_Med_basins(ARGO)
           VALCLIM  = float(df_clim.loc[df_clim.index==NAME_BASIN].iloc[:,COLINDEX])
           if tmp_meta.TREND_TIME_SERIES.isnull().values.any():
              OFFSET  = np.float(tmp.VAR) - VALCLIM
           else:
              Corrrected_val = np.float(tmp.VAR) - np.float(tmp_meta.TREND_TIME_SERIES)
              OFFSET  = Corrrected_val - VALCLIM
           STDCLIM   = float(df_cstd.loc[df_cstd.index==NAME_BASIN].iloc[:,COLINDEX])
           STDCLIM_2 = 2*STDCLIM
           df_report.loc[ df_report.WMO == WMO , 'OFFSET'] = OFFSET
           df_report.loc[ df_report.WMO == WMO , 'basin'] = NAME_BASIN
           df_report.loc[ df_report.WMO == WMO , 'EMODNET'] = VALCLIM
           df_report.loc[ df_report.WMO == WMO , 'STdev*2'] =  STDCLIM_2
           if abs(OFFSET) >= STDCLIM_2:
              print ('rejected data --> black list')
              df_report.loc[ (df_report.WMO == WMO) & (df_report.Depth== DEPTH), 'Black_list'] = 'True'
              OUTPATH_NAME =  OUT_META+ "Blacklist_wmo.csv"
              indexlenght  = 1
              columns_list = ['WMO', 'DATE_DAY' , 'OFFSET' , 'STDCLIM_2']
              values_list  = [int(WMO), DATE_DAY, OFFSET , STDCLIM_2]
              commons_.save_report(OUTPATH_NAME,indexlenght,columns_list, values_list  )
           else:
              print('accepted data -->' +str(WMO) )
              Mask = [Pres<=DEPTH]
              Prof_Coriolis = Profile.copy()

              z_depth = np.array([Pres[0], DEPTH])
              x_var   = np.array([0, np.float(tmp_meta.TREND_TIME_SERIES)])
              if np.isnan(np.sum(x_var)): 
                 pass
              else:
                 fittval = np.polyfit(z_depth, x_var, 1)
                 polyval = np.polyval(fittval, Pres[Mask]  )
                 Profile[0:len(polyval)] = Profile[Mask]   -   polyval
                 Mask_deep = [Pres>DEPTH]
                 Profile[len(polyval):]  = Profile[Mask_deep] - polyval[-1]  
              outname  = superfloat_oxy_function.get_outfile(NPROF, OUTDIR)
              df_report.loc[ df_report.WMO == WMO , 'filename'] =  outname
              Condition_to_write=True
              if Condition_to_write:
                 SAVEDIR = OUTDIR + WMO 
                 if not os.path.exists(SAVEDIR):
                    os.makedirs(SAVEDIR)
                 import superfloat_generator
                 writing_mode=superfloat_generator.writing_mode(outname)
                 condition_to_write = superfloat_generator.exist_valid_variable('DOXY',outname)
                 metadata = NPROF._my_float
                 superfloat_oxy_function.dump_oxygen_file(outname, NPROF , Pres, Profile, Qc, metadata, mode=writing_mode )
                 
Saving_report=True
if Saving_report:
   df_report.to_csv(OUT_META+DATE_DAY+'_report.csv')



import pandas as pd
import matplotlib.pyplot as plt
from basins.region import Region, Rectangle
import basins.V2 as basV2
import numpy as np
import os
import datetime as datetime

def Profile_plot(fig,  x , z , y , label_1, label_2):
  """ lineplot comparing 2 lines x , y
      having same z axis eg. x,y = profiles 1 and 2 and z=depth
  """
  plt.plot(x , z , 'r' , label = label_1)
  plt.plot(y , z , 'b' , label = label_2)
  plt.grid()
  plt.legend()
  return(fig)

def timeseries_plot(fig,  x , y , z , label_1, label_2):
  """ lineplot comparing 2 lines y1 , y2
      having same x axis eg. timeseries
  """
  plt.plot(x , y , 'r' , label = label_1)
  plt.plot(x , z , 'b' , label = label_2)
  plt.grid()
  plt.legend(loc='upper left')
  return(fig)

def col_to_dt(df,name_col):
    """ from object to datetime --> name_col
    """
    df['date'] = pd.to_datetime(df[name_col]).dt.date
    df['date'] = pd.to_datetime(df['date']).dt.date
    df['date'] = pd.to_datetime(df['date'])
    return(df)

def drop_rows_df(df, LIST_ROWS, list_name):
    """ remove specific rows from dataframe 
        by using the list 
    """
    for ITEM in LIST_ROWS:
        df = df.drop(index=ITEM)
    df.index=list_name
    return(df)

def cross_Med_basins(RECTANGLE):
   if RECTANGLE.cross(basV2.eas3):
      LIST_REGION=['lev1','lev2','lev3','lev4','aeg']
      if RECTANGLE.cross(basV2.lev1):
         return(basV2.lev1.name, basV2.lev1.borders, )
      elif RECTANGLE.cross(basV2.lev2):
         return(basV2.lev2.name, basV2.lev2.borders)
      elif RECTANGLE.cross(basV2.lev3):
         return(basV2.lev3.name, basV2.lev3.borders)
      elif RECTANGLE.cross(basV2.lev4):
         return(basV2.lev4.name, basV2.lev4.borders)
      elif RECTANGLE.cross(basV2.aeg):
         return(basV2.aeg.name, basV2.aeg.borders)
   elif  RECTANGLE.cross(basV2.wes3):
      LIST_REGION=['alb','nwm','tyr1','tyr2','swm1','swm2']
      if RECTANGLE.cross(basV2.alb):
         return(basV2.alb.name, basV2.alb.borders)
      elif RECTANGLE.cross(basV2.nwm):
         return(basV2.nwm.name, basV2.nwm.borders)
      elif RECTANGLE.cross(basV2.tyr1):
         return(basV2.tyr1.name, basV2.tyr1.borders)
      elif RECTANGLE.cross(basV2.tyr2):
         return(basV2.tyr2.name, basV2.tyr2.borders)
      elif RECTANGLE.cross(basV2.swm1):
         return(basV2.swm1.name, basV2.swm1.borders)
      elif RECTANGLE.cross(basV2.swm2):
         return(basV2.swm2.name, basV2.swm2.borders)
   elif RECTANGLE.cross(basV2.mid3):
      LIST_REGION=['adr1','adr2','ion1','ion2','ion3']
      if RECTANGLE.cross(basV2.adr1):
         return(basV2.adr1.name, basV2.adr1.borders)
      elif RECTANGLE.cross(basV2.adr2):
         return(basV2.adr2.name, basV2.adr2.borders)
      elif RECTANGLE.cross(basV2.ion1):
         return(basV2.ion1.name, basV2.ion1.borders)
      elif RECTANGLE.cross(basV2.ion2):
         return(basV2.ion2.name, basV2.ion2.borders)
      elif RECTANGLE.cross(basV2.ion3):
         return(basV2.ion3.name, basV2.ion3.borders)

def time_serie_plot(LIST_DEPTH, WMO_LIST, df , df_report):
   for DEPTH in LIST_DEPTH:
     #df['Corrected'] = np.nan
     tmp_data = df.loc[df.Depth==DEPTH]
     tmp_meta = df_report[df_report.Depth==DEPTH]
     for WMO in WMO_LIST:
        tmp_data = tmp_data.loc[tmp_data.name==WMO]
        tmp_data['Corrected'] = np.nan
        tmp_meta = df_report[df_report.WMO==WMO]
        if tmp_meta.DRIFT_CODE.iloc[0] ==1:
           LIST_AVG_TREND=[]
           for DDay in range(0,len(tmp_data)):
              tmp_data.time = pd.to_datetime(tmp_data.time)
              tmpDay = tmp_data.time.iloc[DDay]
              ts_dur = tmpDay- tmp_data.time.min()
              ts_day = ts_dur.days
              TREND_TIME_SERIES = tmp_meta['TREND_TIME_SERIES'].iloc[0]
              DAYS   = tmp_meta['DURATION'].iloc[0]
              TMP_SLOPE = (TREND_TIME_SERIES*ts_day)/DAYS #
              # TREND_TIME_SERIES : DAYS = TMP_SLOPE : ts_day
              LIST_AVG_TREND.append(TMP_SLOPE)
           tmp_data['Corrected'] = np.array(tmp_data.loc[tmp_data.name==WMO]['VAR'] - LIST_AVG_TREND)
           tmp_data = col_to_dt(tmp_data,'time')
           fig, axs = plt.subplots(1, figsize=(9, 5.5))
           timeseries_plot(fig, np.array(tmp_data.time) , np.array(tmp_data.VAR), np.array(tmp_data.Corrected) ,  'Coriolis', 'Corrected')
           plt.tick_params(labelrotation=45)
           plt.suptitle('Coriolis ARGO correction WMO: ' + str(WMO) +' at '+str(DEPTH)+ ' with trend =' + str(np.float(np.around(tmp_meta.TREND_per_YEAR.iloc[0],2)) ) + ' mmol/m3/yr')
           plt.savefig(str(WMO) +'_'+str(DEPTH)+ '_TimeSeries' )


def save_report(OUTPATH_NAME, indexlenght, columns_list, values_list):
   """ OUTPATH_NAME EG. 'OUTPUTS/High_time_freq_argo.csv'
   """ 
   if os.path.exists(OUTPATH_NAME):
       df  = pd.read_csv(OUTPATH_NAME,index_col=0)
   else:
       df  = pd.DataFrame(index=np.arange(0, indexlenght), columns= columns_list  )
   dftmp = pd.DataFrame(index=np.arange(0,1), columns= columns_list )
   dftmp = pd.Series(values_list , columns_list) 
   df =  df.append(dftmp , ignore_index=True)
   df.drop_duplicates(inplace=True)
   df = df.sort_values(by="DATE_DAY")
   df.to_csv(OUTPATH_NAME)





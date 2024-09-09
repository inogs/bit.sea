from basins import V2 as basV2
import os
import pandas as pd
import numpy as np
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

def save_report(OUTPATH_NAME, indexlenght, columns_list, values_list):
   """ OUTPATH_NAME EG. 'OUTPUTS/High_time_freq_argo.csv'
   """
   if os.path.exists(OUTPATH_NAME):
       df  = pd.read_csv(OUTPATH_NAME,index_col=0)
   else:
       df  = pd.DataFrame(index=np.arange(0, indexlenght), columns= columns_list  )
   dftmp = pd.DataFrame(index=np.arange(0,1), columns= columns_list )
   dftmp = pd.Series(values_list , columns_list)
   df = pd.concat((df, dftmp),ignore_index=True)
   df.drop_duplicates(inplace=True)
   df = df.sort_values(by="DATE_DAY")
   df.to_csv(OUTPATH_NAME)
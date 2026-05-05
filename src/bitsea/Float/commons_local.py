from bitsea.basins import V2 as basV2
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

def save_report(lock_fd, indexlenght, columns_list, values_list):
   """ lock_fd EG. open('OUTPUTS/High_time_freq_argo.csv', 'a+')
   """
   lock_fd.seek(0, os.SEEK_END)
   file_is_empty = lock_fd.tell() == 0
   if file_is_empty:
      df = pd.DataFrame(index=np.arange(0, indexlenght), columns=columns_list)
   else:
      lock_fd.seek(0)
      df = pd.read_csv(lock_fd, index_col=0)
   df.dropna(how="all", inplace=True)
   #print(f"loaded df: size: {df.size}; shape: {df.shape}")
   dftmp = pd.Series(values_list, columns_list)
   #print(f"tmp series: size: {dftmp.size}; shape: {dftmp.shape}")
   df = pd.concat((df, dftmp.to_frame().T), ignore_index=True)
   #print(f"after concat df: size: {df.size}; shape: {df.shape}")
   df.drop_duplicates(inplace=True)
   #print(f"dropped duplicates df: size: {df.size}; shape: {df.shape}")
   #print("DATAFRAME AFTER CONCAT AND DROPPING DUPLICATES:")
   #print(df)
   df = df.sort_values(by="DATE_DAY")
   #print(f"sorted df: size: {df.size}; shape: {df.shape}")

   lock_fd.seek(0)
   lock_fd.truncate()
   df.to_csv(lock_fd)
   lock_fd.flush()
   os.fsync(lock_fd.fileno())
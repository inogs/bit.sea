import pandas as pd 
import numpy as np
import copy
import sys 

def check_data(LIST, day):
    if not LIST:
       sys.exit('No data in the selected for:  ' + day) # no data in the selected day
    else:
       print (' ...    available argo in date ' + day + '   ' + str(LIST))

def check_lenght_timeseries(df):
    """input:  dataframe 
       output: Boolean True --> time series < year
       Method: it uses pandas datetime to count days of timeseries  

       """
    df1=df.copy()
    df1.dropna(inplace=True)
    df1.sort_values(by=['time'])
    df1['TIME']  =  pd.to_datetime(df1.time)
    df1.date = pd.to_datetime(df1.time)
    min_d = df1.date.min()
    max_d = df1.date.max()
    duration = max_d-min_d
    days = duration.days
    return(days>365)

# check over nans
def nans_check(df, namecol):
    """ input  : dataframe & nameof column selected for the check 
        output : Boolean
        method : fill nans inside the dataset to remove the starting part of time series: create list_index
                 nan are counted in the core of time series df_ix[list_index]
        """
    df[namecol]   = pd.to_numeric(df[namecol])
    list_index    = list((df.interpolate()).dropna().index)
    return(df.ix[list_index].VAR.isna().sum()  < df.shape[0]/2)

def lenght_timeseries(df, namecol):
    df1=df.copy()
    df1.dropna(inplace=True)
    df1.sort_values(by=[namecol])
    df1[namecol]  =  pd.to_datetime(df1[namecol])
    min_d = df1[namecol].min()
    max_d = df1[namecol].max()
    duration = max_d-min_d
    days = duration.days
    return(days ,min_d , max_d)

def rep(LIST_DEPTH, df, WMO_LIST):
    for DEPTH in LIST_DEPTH:
        tmp = df.loc[df.Depth == DEPTH]
        df_rep = rep0(WMO_LIST, tmp)
        result=df_rep.append(df_rep)
        return(result) 

def rep0(WMO_LIST, df):
    Cday , Chalf = 0, 0
    Cdel , Cadj , Cand = 0 , 0 , 0
    LIST_COL  =['tot_argo','argo_less_1yrs','argo_half_nan','Adj','Del','Adj_Del']
    SELECTED_WMO = []
    for WMO in WMO_LIST:
      tmp = df.loc[df.name==WMO]
      deep = tmp.copy()
      Bool = check_lenght_timeseries(deep)
      if not Bool:
         Cday+=1
         pass
      else:
         del (deep, Bool)
         deep = tmp.copy()
         Bool = nans_check(deep, 'Depth')
         if not Bool:
            Chalf+=1
            pass
         else:
            SELECTED_WMO.append(WMO)
            if len(tmp.Type.drop_duplicates()) ==1 and tmp.Type.drop_duplicates().values=='D':
                Cdel+=1
            elif  len(tmp.Type.drop_duplicates()) ==1 and tmp.Type.drop_duplicates().values=='A':
                Cadj+=1
            else:
                len(tmp.Type.drop_duplicates()) ==2
                Cand+=1
    list_val  = [len(WMO_LIST), Cday, Chalf ,Cadj,Cdel,Cand ]
    df_report =  pd.DataFrame(index=np.arange(0,1), columns= LIST_COL)
    df_report = pd.Series(list_val, df_report.columns)
    return(df_report)

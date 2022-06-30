import pandas as pd
df =pd.read_csv('OUTPUTS/Report_20121231_20220103.csv',index_col=0)
df1 = df[df.basin=='nwm']

#rej = df1[df1.Black_list==True]
#acep = df1[df1.Black_list==False]
#rej.DRIFT_CODE[rej.DRIFT_CODE==1]

import numpy as np
df1['new_col'] = np.nan
df1.loc[df1['Black_list'] == True,  'new_col'] =0
df1.loc[(df1['Black_list']==False) & (df1['DRIFT_CODE']==-1) ,'new_col'] = -1
df1.loc[(df1['Black_list']==False) & (df1['DRIFT_CODE']==1) ,'new_col'] = 1
#plot = df1.plot.pie(y='new_col', figsize=(5, 5))
df1['sum'] = 1
import matplotlib as mpl
import matplotlib
import matplotlib.pyplot as plt
#df1.groupby(['new_col']).sum().plot.pie(y='sum', figsize=(5, 5))


import numpy as np

df1 = col_to_dt(df1, 'max_date')
df1['Year']=df1.date.dt.year
df1['month']=df1.date.dt.month





import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression, TheilSenRegressor
from sklearn.linear_model import RANSACRegressor
from sklearn.linear_model import HuberRegressor, Ridge
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import time

def Trend_RANSAC_ThSEN(x, y):
    estimators = [
                 ('Theil-Sen', TheilSenRegressor(random_state=42)),
                 ('RANSAC', RANSACRegressor(random_state=42)),
                 ]
    X = x[:, np.newaxis]
    line_x = np.array([x.min(), x.max()])
    LIST_APP_TRENDS=[]
    for name, estimator in estimators:
        t0 = time.time()
        estimator.fit(X, y)
        elapsed_time = time.time() - t0
        y_pred = estimator.predict(line_x.reshape(2, 1))
        yP     = np.round(y_pred, 2)
        name_est = name
        LIST_APP_TRENDS.append([name_est, y_pred[1]-y_pred[0]])
    return(LIST_APP_TRENDS)

def sign_analysis(ARRAY):
    """ input np.array
        output Boolean
        action: if all item in np.array have same sign True else False
    """
    i = 0
    while i < len(ARRAY)-1:
          if  np.sign(ARRAY[i]) ==  np.sign(ARRAY[i+1]):
              i+=1
              if i == len(ARRAY)-1:
                 return (True)
          else:
                 return(False)


def compute_trend(df):
    AVG_TREND = np.mean( np.array([ np.array(df['Theil-Sen']), np.array(df['RANSAC']) ]))
    YR_TREND  = (AVG_TREND*365)/df.DURATION.max()
    return(AVG_TREND, YR_TREND)



def trend_conditions(WMO,  days, Bool, DEPTH, min_d, max_d , TI_3, tmp):
    """  compute the trend if time-nans and vals for the selected day
         are satisfied. 
         CASE_* is a Boolean 
    """
    CASE_0 = (days >= 365) and (Bool) and (DEPTH == 600) and (  (max_d.date() - TI_3.end_time.date()).days ==0 )
    CASE_1 = (days >= 365) and (Bool) and (DEPTH == 800) and (abs((max_d.date() - TI_3.end_time.date()).days) <=15)
    if CASE_0 or CASE_1:
       y = np.array(tmp.VAR)
       x = np.array(tmp.index)
       LIST_APP_TRENDS = Trend_RANSAC_ThSEN(x , y)
       lst   = [WMO ,DEPTH,  days , min_d , max_d, LIST_APP_TRENDS[0][1], LIST_APP_TRENDS[1][1] ]
    else: 
       lst   = [WMO ,DEPTH,  days , min_d , max_d, np.nan , np.nan ]
    

    return (lst)

def drift_coding(WMO,Bool ,df0,df1):

    AVG_TREND, YR_TREND = compute_trend(df0)
    CASE_0 = (Bool) and (abs(YR_TREND) >= 1)
    CASE_1 = (Bool) and (abs(YR_TREND) <  1)
    CASE_2 = (Bool)
    if CASE_0:
        df1['TREND_TIME_SERIES'].loc[df1.WMO==WMO] = AVG_TREND
        df1['TREND_per_YEAR'].loc[df1.WMO==WMO] = YR_TREND
        df1['DRIFT_CODE'].loc[df1.WMO==WMO] = 1
    elif CASE_1:
        df1['DRIFT_CODE'].loc[df1.WMO==WMO] = 0
        df1['TREND_per_YEAR'].loc[df1.WMO==WMO] = YR_TREND
    else:
        df1['DRIFT_CODE'].loc[df1.WMO==WMO] =-1
    return(df1)


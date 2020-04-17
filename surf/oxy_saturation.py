import numpy as np
import math

def oxy_sat(p):
    '''
    Arguments:
    * p * a profile object
    Calculate the oxygen at saturation with the formula from
    Garcia and Gordon, 1992 L&O
    using Temperature and Salinity "at surface" misurated by
    argo float. Instead of the "real surface", for T and S we adopt
    the mean value of the first 5m.
    
    Returns:
    * O2o * a concentration of oxygen in mmol/m3
    '''

    
    PresT,Temp,QcT=p.read("TEMP")
    PresS,Sal,QcS=p.read("PSAL")

    temp=np.mean(Temp[PresT<=5]) 
    sal =np.mean(Sal[PresS<=5])

    ts = math.log( ( 298.15 - temp) / ( temp + 273.15 ) )
    ts2 = ts * ts 

    lnk = 2.00856 + 3.22400 *ts + 3.99063 *ts2 + 4.80299 *ts *ts2 + 9.78188 * 10**(-01) * ts2 *ts2 + 1.71069*ts2 *ts2 *ts +  sal  * (-6.24097 * 10**(-03) -6.93498 * 10**(-03) * ts -6.90358 * 10**(-03) * ts2 -4.29155 * 10**(-03) *ts *ts2) -3.11680 * 10**(-7)  *sal *sal
    o2sat_ml= np.exp(lnk)  # oxygen saturation in ml / l

    o2sat_mmol=  o2sat_ml * 1000 / 22.391  # Oxygen saturation in (mmol/m3)  % ICES conversion volume of 1 mole O2 at STP

    return o2sat_mmol

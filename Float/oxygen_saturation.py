import numpy as np

depth_threshold = 5.0 #m

def oxy_sat(p):
    '''
    Calculate the oxygen at saturation with the formula from
    Garcia and Gordon, 1992 L&O
    using Temperature and Salinity "at surface" misurated by
    argo float. Instead of the "real surface", for T and S we adopt
    the mean value of the first 5m.
    As input, it requires the profile object.
    Returns:
    * O2o * a concentration of oxygen in mmol/m3
            nan if there are no TEMP and PSAL values with pressure less than 5m.
    '''
    
    if p.has_adjusted:
        PresT, Temp, QcT = p.read('TEMP', read_adjusted=False)#p._my_float.adjusted('TEMP'))
        PresS, Sali, QcS = p.read('PSAL', read_adjusted=False)#p._my_float.adjusted('PSAL'))
    else:
        PresT, Temp, QcT = p.read('TEMP')
        PresS, Sali, QcS = p.read('PSAL')

    ii = PresT<=depth_threshold
    if ii.sum() == 0:
        print "No TEMP values with pressure less than %f m " %(depth_threshold)
        return np.nan

    temp=np.mean(Temp[ii]) 
    sal =np.mean(Sali[ii])

    ts = np.log( ( 298.15 - temp) / ( temp + 273.15 ) )
    ts2 = ts * ts 

    lnk = 2.00856 + 3.22400 *ts + 3.99063 *ts2 + 4.80299 *ts *ts2 + 9.78188 * 10**(-01) * ts2 *ts2 + 1.71069*ts2 *ts2 *ts +  sal  * (-6.24097 * 10**(-03) -6.93498 * 10**(-03) * ts -6.90358 * 10**(-03) * ts2 -4.29155 * 10**(-03) *ts *ts2) -3.11680 * 10**(-7)  *sal *sal
    o2sat_ml= np.exp(lnk)  # oxygen saturation in ml / l

    o2sat_mmol=  o2sat_ml * 1000 / 22.391  # Oxygen saturation in (mmol/m3)  % ICES conversion volume of 1 mole O2 at STP

    return o2sat_mmol



def oxy_check(Pres,Prof,p):
    '''
    Check float oxygen quality comparing surface value with oxygen saturation
    Argument:
    * Pres * numpy 1d array, pressure profile
    * Prof * Oxygen values profile
    * p * float profile object (needed for T and S inn the Oxy Sat calculation)

    Return:
    * CR *  logical value, True if the check is passed
            False when oxy_Sat returns nan
            False when min(Pres) < 5m
    
    
    '''


    Osat = oxy_sat(p)
    if np.isnan(Osat): return False

    bad=(np.isnan(Prof))
    Prof_good = Prof[~bad]
    Pres_good = Pres[~bad]

    ii = Pres_good<=depth_threshold
    if ii.sum() == 0 :
        print "No DOXY values with pressure less than %f m " %(depth_threshold)
        return False
# Define O2o surface as the mean values in the first 5m from the sea surface:
    Osurf=np.mean(Prof_good[ii])

    mydiff = np.abs(Osurf-Osat)
    
    if (mydiff < 10):
       return True
    else:
       return False


if __name__ == "__main__":

    from commons.time_interval import TimeInterval
    from instruments import bio_float
    from basins import V2 as OGS


    DATESTART = "20190101"
    DATE__END = "20191231"

    T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d')
    var="DOXY"

    Profilelist=bio_float.FloatSelector(var,T_INT,OGS.med)
    p=Profilelist[100]

    Pres, Prof, Qc= p.read(var,read_adjusted=True)
    CR=oxy_check(Pres,Prof,p)

    print "DOXY check is " + str(CR)

    if not CR:
               Prof=np.nan*np.ones_like(Prof)
               Qc=4*np.ones_like(Qc)
               print "Exclude WMO " + np.str(p._my_float.wmo) + " - " + p._my_float.time.strftime("%Y%m%d")


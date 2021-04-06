import numpy as np
from validation.deliverables.float_OXY_saturation import oxy_sat

def oxy_check(Pres,Prof,p):
    '''
    Check float oxygen quality comparing surface value with oxygen saturation
    Argument:
    * Pres * numpy 1d array, pressure profile
    * Prof * Oxygen values profile
    * p * float profile object (needed for T and S inn the Oxy Sat calculation)

    Return:
    * CR * Check Result wich is 0 ==> check passed
                                1 ==> check failed
    '''


    Osat = oxy_sat(p)

    bad=(np.isnan(Prof))
    Prof_good = Prof[~bad]
    Pres_good = Pres[~bad]
 
# Define O2o surface as the mean values in the first 5m from the sea surface:
    Osurf=np.mean(Prof_good[Pres_good<=5])

    mydiff = np.abs(Osurf-Osat)
    
    if (mydiff < 10):
       CR=0 # CHECK RESUT = good
    else:
       CR=1 # CHECK RESULT = bad
    # print mydiff

    return CR

if __name__ == "__main__":

    from commons.time_interval import TimeInterval
    from commons.Timelist import TimeList
    import matplotlib.pyplot as pl
    import matplotlib.dates as mdates
    import numpy.ma as ma
    from instruments import superfloat as bio_float
    from basins import V2 as OGS

    from instruments.check_float_oxygen import oxy_check

    DATESTART = "20190101"
    DATE__END = "20191231"

    T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d')
    var="DOXY"

    Profilelist=bio_float.FloatSelector(var,T_INT,OGS.med)
    p=Profilelist[100]

    Pres, Prof, Qc= p.read(var)
    CR=oxy_check(Pres,Prof,p)

    print "DOXY check is " + str(CR)

    if (CR==1):
               Prof=np.nan*np.ones_like(Prof)
               Qc=4*np.ones_like(Qc)
               print "Exclude WMO " + np.str(p._my_float.wmo) + " - " + p._my_float.time.strftime("%Y%m%d")


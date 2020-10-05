import numpy as np
from datetime import datetime
from scipy import interpolate
import os,sys
from commons.utils import addsep

from instruments import bio_float as bio_float
from instruments.var_conversions import FLOATVARS
from datetime import datetime
from Float import canyon_b_N3n
from Float.canyon_b_N3n import *
#from Float.canyon_b_N3n import get_nitrate


def canyon_nitrate_correction(p):
    
    "from the profile object, get the correction using canyon routine"

    Tp, T, Tqc=p.read('TEMP', False)
    Sp, S, Sqc = p.read('PSAL', False)
    Np, N, Nqc = p.read(FLOATVARS['N3n'],True)
    OXp, OX, OXqc = p.read(FLOATVARS['O2o'],True)

    # Check value at 900m depth

    t_lev = 900  # target level
    if (Np[-1] >= t_lev):
        ii = Np>=t_lev
    elif (Np[-1]<t_lev) & (Np[-1]>100):
        t_lev = Np[-1]
        ii = Np>=t_lev
#    else: continue

    print t_lev   
    iiOX = OXp>=t_lev
    iiT = Tp>=t_lev


    p900=Np[ii][0]
    N900=N[ii][0]
    O900=OX[iiOX][0]
    T900=T[iiT][0]
    S900=S[iiT][0]    
    prof_time=p.time.strftime("%Y%m%d-%H:%M:%S")
    prof_lat=p.lat
    prof_lon=p.lon


    d=datetime.strptime(prof_time,"%Y%m%d-%H:%M:%S")
    lat=prof_lat
    lon=prof_lon
    pres=p900
    temp=T900
    psal=S900
    doxy=O900

    nit=get_nitrate(d, lat, lon, pres, temp, psal, doxy)

    shift=N900-nit

    New_N=N-shift


    ii = (New_N < 0) & (Np < 200)
    New_N[ii] = 0.05
    ii = New_N > 0
    Np = Np[ii]
    New_N = New_N[ii]
    Nqc   =   Nqc[ii]

    return New_N, Np, Nqc



if __name__ == "__main__":
    gtime="20141209-08:45:00"
    lat=17.6; lon=-24.3; pres=float(180); temp=16; psal=36.1; doxy=104
    d=datetime.strptime(gtime,"%Y%m%d-%H:%M:%S")
    for i in range(10): print get_nitrate(d, lat, lon, pres, temp, psal, doxy)

\

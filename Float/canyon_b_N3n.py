import numpy as np
from datetime import datetime
import os,sys
from commons.utils import addsep
from mhelpers.linear_shift import linear_shift

Training_dir="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/CANYON-B/"
basedir=addsep(os.getenv("CANYONB_TRAINING_DIR", Training_dir))
if not os.path.exists(basedir):
    print basedir
    raise ValueError("Environment variable CANYONB_TRAINING_DIR must be defined")

inwgts=np.loadtxt(basedir + 'wgts_NO3.txt')


def get_nitrate(timeobj,lat,lon, pres,temp, psal, doxy):
    '''
    Calculates the nitrate value of canyon_med
    Arguments:
    * timeobj * a datetime object
    * lat     * scalar value
    * lon     * idem
    * pres    * idem
    * psal    * idem
    * doxy    * idem
    '''

    lon = float(lon)
    lat = float(lat)
    pres = float(pres)
    temp = float(temp)
    psal = float(psal)
    doxy = float(doxy)
    doy=np.int(timeobj.strftime('%j'))*360./365
    
    if lon>180: lon=lon-360
    data=[lat/90,  np.abs(1-np.mod(lon-110,360)/180.), np.abs(1-np.mod(lon-20,360)/180.), temp,psal,doxy,pres/2e4+1./(1+np.exp(-pres/300.))**3 ]
    data = np.array(data)
    no = 1
    nol= 1


    noparsets=inwgts.shape[1]-1
    ni = len(data)
    mw=inwgts[0:ni+1,-1]
    sw=inwgts[ni+1:2*ni+2,-1]
    data_N=(data-mw[:ni])/sw[:ni]
    wgts=inwgts[3,0:noparsets]
    betaciw=inwgts[2*ni+2:,noparsets]
    ii=np.isnan(betaciw)
    betaciw=betaciw[~ii]


    cval   = np.ones((nol,noparsets), np.float32)*np.nan
    cvalcy = np.ones((  1,noparsets), np.float32)*np.nan
    inval  = np.ones((nol,ni,noparsets), np.float32)*np.nan

    for l in range(noparsets):
        nlayerflag=1+np.bool(inwgts[1,l])
        nl1=int(inwgts[0,l])
        nl2=int(inwgts[1,l])
        beta=inwgts[2,l]
        pos_1=nl1*ni+4
        pos_2 = pos_1 + nl1
        pos_3 = pos_2 + nl2*nl1
        pos_4 = pos_3 + nl2
        w1 = inwgts[4:pos_1,l].reshape(ni,nl1).T
        b1 = inwgts[pos_1:pos_2,l]
        w2 = inwgts[pos_2:pos_3,l].reshape(nl1,nl2).T
        b2 = inwgts[pos_3:pos_4,l]

        if nlayerflag ==2 :
            pos_5 = pos_4 + no*nl2
            pos_6 = pos_5 + no
            w3 = inwgts[pos_4:pos_5,l].reshape(nl2,no).T
            b3 = inwgts[pos_5:pos_6,l]
        a=np.dot(data_N.T,w1.T)+b1
        if nlayerflag ==1:
            y = np.dot(np.tanh(a.T), w2.T) + b2
        if nlayerflag ==2:
            b = np.dot(np.tanh(a.T), w2.T) + b2
            y = np.dot(np.tanh(b.T), w3.T) + b3

        cval[:,l]=y
        cvalcy[:,l] = 1./beta
        x1 = w1 * (1-np.tanh(a)**2).reshape(nl1, 1)
        if nlayerflag == 1:
            inx = np.dot(w2,x1)
        if nlayerflag == 2:
            x2 = w2 * (1-np.tanh(b)**2).reshape(nl2,1)
            inx = np.dot(w3,x2).dot(x1)
        inval[:,:,l] = inx


    # Denormalization of the network output
    cval=cval*sw[ni]+mw[ni] # variable
    cvalcy=cvalcy*sw[ni]**2 # 'noise' variance
    # add committee of all params as evidence-weighted mean
    V1=wgts.sum()
    V2=(wgts**2).sum()
    out_value = (wgts*cval).sum()/V1 # weighted mean
    return out_value

    cvalcu = (wgts*(cval-out_value)**2).sum()/(V1 - V2/V1) # CU variance


def canyon_nitrate_correction(p, Np, N, Nqc, OXp, OX):

    '''Get the correction using canyon routine

    Arguments:
    * p  *  an object profile,
    * Np *  numpy 1d array, nitrate pressure
    * N  *  numpy 1d array, nitrate values
    * Nqc*  numpy 1d array, nitrate Qc
    * OXp*  numpy 1d array, oxygen pressure
    * OXp*  numpy 1d array, oxygen values

    Returns:
    * Pres  * numpy 1d array, corrected pressure
    * Value *
    * Qc    *
    * t_lev * in meters, target level of the correction
    * nit   * nitrate value at the target level
    '''



    Tp, T, _=   p.read('TEMP', False)
    Sp, S, _  = p.read('PSAL', False)
#    Np, N, Nqc = p.read('NITRATE',True)
#    OXp, OX, _ = p.read('DOXY',True)

    # Check value at 900m depth

    t_lev = 900  # target level
    if (Np[-1] >= t_lev):
        ii = Np>=t_lev
    elif (Np[-1]<t_lev) & (Np[-1]>100):
        t_lev = Np[-1]
        ii = Np>=t_lev
#    else: continue


    iiOX = OXp>=t_lev
    iiT = Tp>=t_lev

    p900=Np[ii][0]
    N900=N[ii][0]
    if iiOX.sum() > 0:
        O900=OX[iiOX][0]
    else:
        O900 = min(t_lev,OXp.max())

    if iiT.sum() > 0:
        T900=T[iiT][0]
        S900=S[iiT][0]
    else:
        T900 = T[-1]
        S900 = S[-1]

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

    # get the nitrate value at "900m" with CanyonB:
    nit=get_nitrate(d, lat, lon, pres, temp, psal, doxy)

    shift=N900-nit
  
    depth_bot=min(600,t_lev)
    New_profile  = linear_shift(N,Np,shift,p_bot=depth_bot)
    Nqc[:] =8

    return Np, New_profile, Nqc, t_lev, nit, shift

if __name__ == "__main__":
    timeobj=datetime(2014,12,9,8,45)
    lat=17.6;lon=-24.3;pres=180.;temp=16;psal=36.1;doxy=104 # test values
    print get_nitrate(timeobj, lat, lon, pres, temp, psal, doxy)
    
    from commons.time_interval import TimeInterval
    from commons.Timelist import TimeList
    import matplotlib.pyplot as pl
    import matplotlib.dates as mdates
    import numpy.ma as ma
    from basins import V2 as OGS
    from instruments import bio_float as bio_float

    from commons.layer import Layer
    from static.climatology import get_climatology

    DATESTART = "20190101"
    DATE__END = "20191231"

    T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d')
    var="NITRATE"

    Profilelist=bio_float.FloatSelector(var,T_INT,OGS.med)
    p=Profilelist[30]

    Pres, Prof, Qc= p.read(var,True)
    print p.available_params
    print p.ID()

    DOXYp, DOXY, QcO= p.read("DOXY",True) 

    print len(DOXY)
    if ((len(Prof)>5) & (len(DOXY)>5)):
    # PERFORM THE CORRECTION:
       Np, N, Qc, t_lev, nit = canyon_nitrate_correction(p, Pres, Prof, Qc, DOXYp, DOXY)

       fig , ax = pl.subplots()
#    ax.plot(Prof,Pres,'b', N, Np,'r')
       ax.plot(Prof,Pres,'b',label="Original")
       ax.plot(N, Np,'r',label="Corrected [Canyon]")

# Check with CLIMATOLOGY:

       SUBLIST = OGS.Pred.basin_list
       PresDOWN=np.array([0,25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500])
       LayerList=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]
       N3n_clim, N3n_std = get_climatology('N3n', SUBLIST, LayerList)
       z_clim = np.array([-(l.bottom+l.top)/2  for l in LayerList])

       for iSub,sub in enumerate(SUBLIST):
        if sub.is_inside(p.lon,p.lat):
           print N3n_clim[iSub,:]
           print sub.name
           sub_name=sub.name
           ax.plot(N3n_clim[iSub,:],-1*z_clim,'g',label="Clim")



       ax.invert_yaxis()
       ax.set_ylabel("depth $[m]$",color = 'k', fontsize=10)
       fig_title="Float ID " + np.str(p._my_float.wmo) + " NITRATE CANYON- " + p._my_float.time.strftime("%Y%m%d")
       ax.set_title(fig_title)
       ax.legend()
#    fig.show()
       fig_name='FLOAT_Canyon_' + np.str(p._my_float.wmo) + '_' + p._my_float.time.strftime("%Y%m%d")  + '.png'
       fig.savefig(fig_name)


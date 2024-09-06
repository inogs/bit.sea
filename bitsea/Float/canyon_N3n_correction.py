import numpy as np
from datetime import datetime
from scipy import interpolate
import os,sys
from commons.utils import addsep

from instruments import bio_float as bio_float
from instruments.var_conversions import FLOATVARS
from datetime import datetime
#from Float.canyonb_N3n import get_nitrate


Training_dir="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/CANYON_B/CODES/CANYON_Training/"
basedir=addsep(os.getenv("CANYONB_TRAINING_DIR", Training_dir))
if not os.path.exists(basedir):
    print basedir
    raise ValueError("Environment variable CANYONB_TRAINING_DIR must be defined")

presgrid=np.loadtxt(basedir + 'CY_doy_pres_limit.csv',delimiter="\t")
poids=np.loadtxt(basedir + 'Fichier_poids_NO3_hidden_ascii20_17.sn')
poids=poids[:,2]
# Mean and standard deviation of the training dataset
# These values are used to normalize the inputs parameters
Moy  =np.loadtxt(basedir + 'moy_NO3.dat')
Ecart=np.loadtxt(basedir + 'std_NO3.dat')




def get_nitrate(timeobj,lat,lon, pres,temp, psal, doxy):
    '''
    Calculates the nitrate value of canyon_b
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
    x,y= np.meshgrid(presgrid[1:,0],presgrid[0,1:])
    f = interpolate.interp2d(presgrid[1:,0],presgrid[0,1:],presgrid[1:,1:].T)
    prespivot = f(lon,lat)

    fsigmoid=1./(1 + np.exp((pres-prespivot)/50))
    lonrad=lon*np.pi/180
    doyrad=doy*np.pi/180

    data=[lat/90, np.sin(lonrad), np.cos(lonrad), np.sin(doyrad)*fsigmoid, np.cos(doyrad)*fsigmoid, temp, psal, doxy, pres/2e4 + 1./((1+np.exp(-pres/300))**3)];
    data = np.array(data)




    ne=9    # Number of inputs
    nc1=20  # Number of neurons of the first hidden layer
    nc2=17  # Number of neurons of the second hidden layer
    ns=1    # Number of outputs


    # WEIGHT AND BIAS PARAMETERS
    # weight and bias from the input layer to the first hidden layer w1, b1
    b1=poids[0:nc1]
    start=nc1+nc2+ns
    end=start+ne*nc1
    w1 = np.reshape(poids[start:end],(ne,nc1)).T

    #% weight and bias from the first hidden layer to the second hidden layer
    b2=poids[nc1:nc1+nc2]
    start=nc1+nc2+ns+ne*nc1
    end = start + nc1*nc2
    w2=np.reshape(poids[start:end],(nc1,nc2)).T
    #% weight and bias from the second hidden layer to the output layer w3, b3
    start=nc1+nc2
    b3=poids[start:start+ns]
    start=nc1+nc2+ns+ne*nc1+nc1*nc2
    end = start + nc2*ns
    w3=np.reshape(poids[start:end],(nc2,ns)).T



    # NORMALISATION OF THE INPUT PARAMETERS
    rx=1 #[rx,~]=size(data);
    data_N=(2./3)*(data-Moy[:ne])/Ecart[:ne]

    # Two hidden layers
    a=1.715905*np.tanh((2./3)*(np.dot(data_N,w1.T)+b1)) # input layer to first hidden layer
    b=1.715905*np.tanh((2./3)*(np.dot(     a,w2.T)+b2)) # first hidden layer to second hidden layer
    y=                         np.dot(     b,w3.T)+b3   # second hidden layer to output layer
    # Y is the normalised output value of the neural network,


    # Denormalisation of the output of the NN for getting the true value
    y_rescaled=1.5*y*Ecart[ne]+Moy[ne]

    # and put into same shape as the input variables
    out=y_rescaled[0]
    return out


#if __name__ == "__main__":
#    gtime="20141209-08:45:00"
#    lat=17.6; lon=-24.3; pres=float(180); temp=16; psal=36.1; doxy=104
#    d=datetime.strptime(gtime,"%Y%m%d-%H:%M:%S")
#    for i in range(10): print get_nitrate(d, lat, lon, pres, temp, psal, doxy)

def canyon_nitrate_correction(p):
    
    "from the profile object, get the correction using canyon routine"

#import numpy as np
#from commons.mask import Mask
#from commons.Timelist import TimeList, TimeInterval
#from instruments import bio_float as bio_float
#from instruments.var_conversions import FLOATVARS
#from basins import V2 as OGS
#from basins.region import Rectangle
#from datetime import datetime
#from Float.canyonb_N3n import get_nitrate

#TI=TimeInterval("20200101","20200201","%Y%m%d")

#ALL_PROFILES = bio_float.FloatSelector(FLOATVARS["N3n"], TI, Rectangle(-6,36,30,46))

#p=ALL_PROFILES[10]

#Np, N, Nqc= p.read(FLOATVARS["N3n"])
    Tp, T, Tqc=p.read('TEMP', False)
    Sp, S, Sqc = p.read('PSAL', False)
    Np, N, Nqc = p.read(FLOATVARS['N3n'],True)
    OXp, OX, OXqc = p.read(FLOATVARS['O2o'],True)

    # Check value at 900m depth
    ii = Np>=900
    iiOX = OXp>=900
    iiT = Tp>=900

    p900=Np[ii][0]
    N900=N[ii][0]
    O900=OX[iiOX][0]
    T900=T[iiT][0]
    S900=S[iiT][0]
    prof_time=p.time.strftime("%Y%m%d-%H:%M:%S")
    prof_lat=p.lat
    prof_lon=p.lon

#gtime="20141209-08:45:00"
#    lat=17.6; lon=-24.3; pres=float(180); temp=16; psal=36.1; doxy=104
#    d=datetime.strptime(gtime,"%Y%m%d-%H:%M:%S")

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
    return New_N



if __name__ == "__main__":
    gtime="20141209-08:45:00"
    lat=17.6; lon=-24.3; pres=float(180); temp=16; psal=36.1; doxy=104
    d=datetime.strptime(gtime,"%Y%m%d-%H:%M:%S")
    for i in range(10): print get_nitrate(d, lat, lon, pres, temp, psal, doxy)

\

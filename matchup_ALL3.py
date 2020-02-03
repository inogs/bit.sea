
from basins import V2 as OGS
from commons.layer import Layer
from instruments import optbio_float_2019
from instruments import var_conversions
import numpy as np
import os,sys
from profiler_ET import *
from profiler_SAT import *
from scipy.optimize import curve_fit
from commons import netcdf4
from instruments.matchup_manager_SAT import Matchup_Manager as Matchup_Manager_SAT

def euphotic(Pres, Ed):
    EUPH_model         = [x for x in range(len(Pres)) if Ed.Model[x]/Ed.Model[0] > 0.37]
    EUPH_float         = [x for x in range(len(Pres)) if Ed.Ref[x]/Ed.Ref[0] > 0.37]
    ind_model          = EUPH_model[-1]
    ind_float          = EUPH_float[-1]
    ind  = ind_float  if ind_float < ind_model else ind_model
    EUPH = EUPH_float if ind_float < ind_model else EUPH_model
    zeu                = Pres[ind]
    Pres_Eu            = Pres[EUPH]
    Ed_Eu_MODEL        = Ed.Model[EUPH]
    Ed_Eu_FLOAT        = Ed.Ref[EUPH]
    return zeu, Pres_Eu, Ed_Eu_MODEL, Ed_Eu_FLOAT


M_model = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
M_sat   = Matchup_Manager_SAT(ALL_PROFILES_sat, TL_sat, BASEDIR_sat)

maskfile    =  '/galileo/home/userexternal/eterzic0/MASKFILE_70/meshmask.nc'
nav_lev = netcdf4.readfile(maskfile,"nav_lev")
jpk = nav_lev.shape[0]
mydepth = nav_lev[0:jpk-1]  

#TI = TimeInterval("20120101", "20171231","%Y%m%d")
TI = TimeInterval("20120101", "20150801","%Y%m%d")

varname     = var_conversions.FLOAT_OPT_BIOPTIMOD['Ed490f'] 


Profilelist    =optbio_float_2019.FloatSelector(varname,TI , OGS.med) 


Kd_Model = 0.1
Kd_Float = 0.1

MODEL = []
FLOAT = []
SAT   = []

for p in Profilelist:
    profile_ID = p.ID()
    print(profile_ID)

    '''
    phase 1. Read BGC-ARGO profiles
    '''
    Pres, Ed_float,  Qc = p.read(varname)   #'IRR_490'
    Lon = p.lon
    Lat = p.lat
    timestr = p.time.strftime("%Y%m%d-%H:%M:%S")
    nLevels = len(Pres)

    if not TI.contains(p.time): continue

    '''
    phase 2. Read model data - in theory you can use this also for float data
    '''
    Ed_matchup_Model = M_model.getMatchups([p], mydepth, 'Ed490f', interpolation_on_Float=True)


    Ed_matchup_SAT   = M_sat.getMatchups( [p],  mydepth, 'KD490' , interpolation_on_Float=False)

    '''
    phase 3. Euphotic depth range calculation
    '''
    zmax, PresEu, Ed_MODEL, Ed_REF = euphotic(Pres, Ed_matchup_Model)

    '''
    phase 4. Calculate Kd
    '''

    Kd_Sat = Ed_matchup_SAT.Model[0]


    if len(PresEu) < 5: continue
    if PresEu[1] - PresEu[0] > 9. : continue
    if PresEu[4] > 15. : continue


    func = lambda z, Kd: np.exp(-Kd*z)

    poptM, pcovM = curve_fit(func, PresEu-PresEu[0], Ed_MODEL/Ed_MODEL[0], p0=Kd_Model)#,bounds=(0.001,0.005))
    poptF, pcovF = curve_fit(func, PresEu-PresEu[0], Ed_REF/Ed_REF[0],     p0=Kd_Float)#,bounds=(0.001,0.005))

    Kd_Model = poptM[0]   ; Kd_Float = poptF[0]

    #MODEL.append(Kd_Model)
    #FLOAT.append(Kd_Float)

    if Kd_Model < 0.01:
        continue

    MODEL.append(Kd_Model)
    FLOAT.append(Kd_Float)
    SAT.append(Kd_Sat)

    print(Kd_Model, Kd_Float, Kd_Sat)
    '''
    phase 4b. Addidionally plot to see if the function works well.
    I'm hereby attaching the script I used in the phase of testing
    '''
    '''
    E0M      = Ed_MODEL[0] ; E0F = Ed_REF[0]
    EdM      = func(PresEu, E0M, Kd_Model)
    EdF      = func(PresEu, E0F, Kd_Float)

    import matplotlib.pyplot as plt
    plt.plot(EdM, -PresEu, 'r')      # modelled (EdF instead of EdM for floats)
    plt.plot(Ed_MODEL, -PresEu, 'g') # data in the euphotic range (Ed_REF instead of Ed_MODEL for floats)
    '''

from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = plt.axes(projection="3d")


z_points = SAT
x_points = MODEL
y_points = FLOAT
ax.scatter3D(x_points, y_points, z_points, c=z_points, cmap='hsv');

ax.set_xlabel('MODEL')
ax.set_ylabel('FLOAT')
ax.set_zlabel('SATELLITE')
plt.show()

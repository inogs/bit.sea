
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


PATH=os.getcwd()   
SIM_NAME = PATH.strip('galileo/home/userexternal/eterzic0/BIOPTIMOD/KD_VAL/INPUT/.../bit.sea/')

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


fid = open(filestat,'wb')

for p in Profilelist[0:20]:
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


    for subbasin in OGS.Pred:
        if subbasin.is_inside(Lon, Lat):
            sub_float = subbasin.name
            print(sub_float)


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

    if Kd_Model < 0.01:
        continue

    
    fid.write("%s %.2f %.2f %.2f %.2f %.2f %.2f \n" % (profile_ID, Kd_Model, Kd_Float, Kd_Sat, Lat, Lon, Sub )


file_dir = PATH + 'STATS_SAT/'
file_out1 =  file_dir + '_' + SIM_NAME.replace('/', '_') + 'test.stat'
fid.close()

    #print(Kd_Model, Kd_Float, Kd_Sat)
     
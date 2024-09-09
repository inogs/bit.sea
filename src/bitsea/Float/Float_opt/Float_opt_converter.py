import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Converts LOV text files of Optical dataset 
    provided by 
    in NetCDF files similar to LOV dataset
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required =False,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V2C/wrkdir/2/POSTPROC/AVE_FREQ_1/TMP/",
                                help = ''' Input directory'''
                                )

    parser.add_argument(   '--outputdir', '-o',
                                type = str,
                                required = True,
                                default = '/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/FLOAT_OPT/',
                                help = 'Operational directory')
    parser.add_argument(   '--lovFloatIndex', '-f',
                                type = str,
                                required = False,
                                default = '/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V5C/FLOAT_LOVBIO/Float_Index.txt',
                                help = 'Operational directory')    

    return parser.parse_args()

args = argument()

import numpy as np
from mydtype import *
from instruments.instrument import ContainerProfile
from datetime import datetime
import scipy.io.netcdf as NC
import os
from commons.utils import addsep

def unique_rows(data, prec=5):
    
    d_r = np.fix(data * 10 ** prec) / 10 ** prec + 0.0
    b = np.ascontiguousarray(d_r).view(np.dtype((np.void, d_r.dtype.itemsize * d_r.shape[1])))
    _, ia,ic = np.unique(b, return_index=True, return_inverse=True)
    return np.unique(b).view(d_r.dtype).reshape(-1, d_r.shape[1]), ia, ic


def get_Profiles(filename, varname, thedtype):
    A=np.loadtxt(filename, dtype=thedtype,skiprows=1)
    lon   =A['long']
    lat   =A['lat']
    year  =A['year']
    month =A['month']
    day   =A['day']
    values=A[varname]
    depth =A['PRES']
    times = year*10000 + month*100 + day
    nValues=len(times)
    M=np.zeros((nValues,3),dtype=np.float64)
    M[:,0] = times
    M[:,1] = lon
    M[:,2] = lat
    uniques,_,ib = unique_rows(M)
    nProfiles=len(uniques)
    Profilelist=[]
    for i in range(nProfiles):
        t = int(uniques[i,0])
        year,yearfrac=divmod(t,10000)
        month,day=divmod(yearfrac,100)
        time = datetime(year,month,day)
        lon  = uniques[i,1].astype(np.float32)
        lat  = uniques[i,2].astype(np.float32)
        inthisprofile = ib==i
        Values=values[inthisprofile]
        Depth = depth[inthisprofile]
        ind = np.nonzero(inthisprofile)[0][0]
        Cruise = A['name_profile'][ind][1:-1]
        LP = ContainerProfile(lon,lat,time,Depth,Values,Cruise)
        Profilelist.append(LP)
    return Profilelist

def Same_Profilelist(p, VARLIST, CHL_PL, ED380_PL, ED412_PL, ED490_PL, PAR_PL, SAL_PL):
    '''
    Returns the a list of profile objects referring to the same lon,lat,time
    The order corresponds to VARLIST
    '''

    Profilelist=[]
    ind =SAL_PL.index(p)
    Profilelist.append(SAL_PL[ind])
    Profilelist.append(p)
    ind = PAR_PL.index(p)
    Profilelist.append(PAR_PL[ind])
    if "CHLA" in VARLIST:
        ind = CHL_PL.index(p)
        Profilelist.append(CHL_PL[ind])
    if "Ed_380" in VARLIST:
        ind = ED380_PL.index(p)
        Profilelist.append(ED380_PL[ind])
    if "Ed_412" in VARLIST:
        ind = ED412_PL.index(p)
        Profilelist.append(ED412_PL[ind])
    if "Ed_490" in VARLIST:
        ind=ED490_PL.index(p)
        Profilelist.append(ED490_PL[ind])
    return Profilelist

INPUTDIR =  addsep(args.inputdir) #"/Users/gbolzon/Downloads/BIO-ARGO-QC-DATA/"
FloatIndexer=args.lovFloatIndex #r"/Users/gbolzon/Downloads/Float_Index.txt"
OUTDIR=addsep(args.outputdir)#"/Users/gbolzon/Downloads/Float_OPT/DATASET/"

if True:
    filename=INPUTDIR + "CHLOROPHYLL_PROFILE_QC_MEDSEA.txt"
    P_chl = get_Profiles(filename, "CHLA_CALIBRATED", CHLOROPHYLL_PROFILE_QC_MEDSEA_type)
    filename=INPUTDIR + "Ed380_PROFILE_TYPE_1_MEDSEA.txt"
    P_Ed380 = get_Profiles(filename, "IRR_380", Ed380_PROFILE_TYPE_1_MEDSEA_type)

    filename=INPUTDIR + "Ed412_PROFILE_TYPE_1_MEDSEA.txt"
    P_Ed412 = get_Profiles(filename, "IRR_412", Ed412_PROFILE_TYPE_1_MEDSEA_type)

    filename=INPUTDIR + "Ed490_PROFILE_TYPE_1_MEDSEA.txt"
    P_Ed490 = get_Profiles(filename, "IRR_490", Ed490_PROFILE_TYPE_1_MEDSEA_type)

    filename=INPUTDIR + "PAR_PROFILE_TYPE_1_MEDSEA.txt" 
    P_PAR = get_Profiles(filename, "IRR_PAR", PAR_PROFILE_TYPE_1_MEDSEA_type)

    filename=INPUTDIR + "SALINITY_PROFILE_MEDSEA.txt"
    P_SAL = get_Profiles(filename, "SAL", SALINITY_PROFILE_MEDSEA_type)
    filename=INPUTDIR + "TEMPERATURE_PROFILE_MEDSEA.txt"
    P_TEM = get_Profiles(filename, "TEMP", TEMPERATURE_PROFILE_MEDSEA_type)



INDEX_FILE=np.loadtxt(FloatIndexer,dtype=FloatIndex_type, delimiter=",",ndmin=1)

UNIQUE_PROFILES=P_TEM[:]

for ip, p in enumerate(UNIQUE_PROFILES):
    for ifile, filename in enumerate(INDEX_FILE['file_name']):
        if filename.find(p.name()) > -1:
            time= INDEX_FILE[ifile]['time'].replace("-","").replace(":","")
            lon = INDEX_FILE[ifile]['lon'].astype(np.float64)
            lat = INDEX_FILE[ifile]['lat'].astype(np.float64)
            break
    else:
        raise ValueError ( "file %s non found" %(p.name()) )

    wmodir=OUTDIR + os.path.dirname(filename)
    outfile=OUTDIR + filename
    os.system("mkdir -p " + wmodir)
    

    VARLIST=['PSAL','TEMP','PAR']
    if p in P_chl  : VARLIST.append('CHLA')   
    if p in P_Ed380: VARLIST.append('Ed_380')
    if p in P_Ed412: VARLIST.append('Ed_412')
    if p in P_Ed490: VARLIST.append('Ed_490')
    
    
    same_profile_objects = Same_Profilelist(p, VARLIST, P_chl,P_Ed380, P_Ed412, P_Ed490, P_PAR, P_SAL)
    
    V=[]
    for ivar, var in enumerate(VARLIST):
        p = same_profile_objects[ivar]
        Pres, Profile, Qc = p.read(var)
        V.append((Pres,Profile,Qc))

# NetCDF generation
    ncOUT = NC.netcdf_file(outfile,"w")
    ncOUT.createDimension("DATETIME",14)
    ncOUT.createDimension("NPROF", 1)
    for ivar, var in enumerate(VARLIST):
        dimname = "n"+var
        dimvalue = len(V[ivar][0])
        ncOUT.createDimension(dimname, dimvalue)
    ncvar=ncOUT.createVariable("REFERENCE_DATE_TIME", 'c', ("DATETIME",))
    ncvar[:]=time
    ncvar=ncOUT.createVariable("JULD", 'd', ("NPROF",))
    ncvar[:]=0.0
    ncvar=ncOUT.createVariable("LONGITUDE", "d", ("NPROF",))
    ncvar[:] = lon
    ncvar=ncOUT.createVariable("LATITUDE", "d", ("NPROF",))
    ncvar[:] = lat
    
    for ivar, var in enumerate(VARLIST):
        ncvar=ncOUT.createVariable("PRES_"+var, 'f', ("n"+var,))
        ncvar[:]=V[ivar][0]
        ncvar=ncOUT.createVariable(var, 'f', ("n"+var,))
        ncvar[:]=V[ivar][1]
    ncOUT.close()


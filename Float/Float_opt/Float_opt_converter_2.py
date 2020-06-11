#modified by eterzic 22.08.2019 from
# bit.sea/Float/Float_opt/Float_opt_converter.py 
# for the BGC-Argo data set from 2012 to 2019

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
    parser.add_argument(   '--FloatIndex', '-f',
                                type = str,
                                required = False,
                                default = '/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V5C/FLOAT_BIO/Float_Index.txt',#'/galileo/home/userexternal/eterzic0/BGC-ARGO-DATA/ORIG/Float_Index.0.txt',
                                help = 'Operational directory')    

    return parser.parse_args()

args = argument()

import numpy as np
from mydtype_2 import *
from instruments.instrument import ContainerProfile
from datetime import datetime
import scipy.io.netcdf as NC
import os
from commons.utils import addsep
from commons.time_interval import TimeInterval

def unique_rows(data, prec=5):
    
    d_r = np.fix(data * 10 ** prec) / 10 ** prec + 0.0
    b = np.ascontiguousarray(d_r).view(np.dtype((np.void, d_r.dtype.itemsize * d_r.shape[1])))
    _, ia,ic = np.unique(b, return_index=True, return_inverse=True)
    return np.unique(b).view(d_r.dtype).reshape(-1, d_r.shape[1]), ia, ic


def get_Profiles(filename, varname="VALUE", dtype=BIOOPTIMOD_type):
    A=np.loadtxt(filename, dtype=dtype,skiprows=1)
    year                            = np.array( [ np.int(A['date'][i][1:5]) for i in range(len(A)) ] )
    month                           = np.array( [ np.int(A['date'][i][6:8]) for i in range(len(A)) ] )
    day                             = np.array( [ np.int(A['date'][i][9:11]) for i in range(len(A)) ] )
    lon   =A['Longitude']
    lat   =A['Latitude']
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
        Cruise = A['Name'][ind][3:10]
        LP = ContainerProfile(lon,lat,time,Depth,Values,Cruise)
        Profilelist.append(LP)
    return Profilelist

#def Same_Profilelist(p, VARLIST, CHL_PL, ED380_PL, ED412_PL, ED490_PL, PAR_PL, SAL_PL):
def Same_Profilelist(p, VARLIST, CHL_PL, ED380_PL, ED412_PL, ED490_PL, PAR_PL):
    '''
    Returns the a list of profile objects referring to the same lon,lat,time
    The order corresponds to VARLIST
    '''

    Profilelist=[]

    if "CHL" in VARLIST:
        ind = CHL_PL.index(p)
        Profilelist.append(CHL_PL[ind])
    if "IRR_380" in VARLIST:
        ind = ED380_PL.index(p)
        Profilelist.append(ED380_PL[ind])
    if "IRR_412" in VARLIST:
        ind = ED412_PL.index(p)
        Profilelist.append(ED412_PL[ind])
    if "IRR_490" in VARLIST:
        ind=ED490_PL.index(p)
        Profilelist.append(ED490_PL[ind])
    if "PAR" in VARLIST:
        ind=PAR_PL.index(p)
        Profilelist.append(PAR_PL[ind])
    return Profilelist

def Same_Profilelist_2020(p, VARLIST, TEMP_PL, SALI_PL, BBP700_PL):
    '''
    Returns the a list of profile objects referring to the same lon,lat,time
    The order corresponds to VARLIST
    '''

    Profilelist=[]

    if "TEMP" in VARLIST:
        ind = TEMP_PL.index(p)
        Profilelist.append(TEMP_PL[ind])
    if "SALI" in VARLIST:
        ind = SALI_PL.index(p)
        Profilelist.append(SALI_PL[ind])
    if "BBP700" in VARLIST:
        ind = BBP700_PL.index(p)
        Profilelist.append(BBP700_PL[ind])
    return Profilelist


INPUTDIR = addsep(args.inputdir)
FloatIndexer=args.FloatIndex
OUTDIR = addsep(args.outputdir)
setup_2019=False
setup_2020=True

if setup_2019:
    filename=INPUTDIR + "QC_CHL_MEDSEA_MAY2019_BIOOPTIMOD.txt"  #"CHL_red.txt"
    P_chl = get_Profiles(filename, "CHL", QC_CHL_MEDSEA_MAY2019_BIOOPTIMOD_type)
    
    filename=INPUTDIR + "QC_380_MEDSEA_MAY2019_BIOOPTIMOD.txt"  #"CHL_red.txt"
    P_Ed380 = get_Profiles(filename, "IRR_380", QC_380_MEDSEA_MAY2019_BIOOPTIMOD_type)
    
    filename=INPUTDIR + "QC_412_MEDSEA_MAY2019_BIOOPTIMOD.txt"  #"CHL_red.txt"
    P_Ed412 = get_Profiles(filename, "IRR_412", QC_412_MEDSEA_MAY2019_BIOOPTIMOD_type)
    
    filename=INPUTDIR + "QC_490_MEDSEA_MAY2019_BIOOPTIMOD.txt"  #"CHL_red.txt"
    P_Ed490 = get_Profiles(filename, "IRR_490", QC_490_MEDSEA_MAY2019_BIOOPTIMOD_type)
    
    filename=INPUTDIR + "QC_PAR_MEDSEA_MAY2019_BIOOPTIMOD.txt"  #"CHL_red.txt"
    P_PAR = get_Profiles(filename, "PAR", QC_PAR_MEDSEA_MAY2019_BIOOPTIMOD_type)
    UNIQUE_PROFILES=P_chl

if setup_2020:
    filename=INPUTDIR + "QC_TEMP_MEDSEA_MAY2020_BIOOPTIMOD.txt"
    P_TEMP = get_Profiles(filename)

    filename=INPUTDIR + "QC_SAL_MEDSEA_MAY2020_BIOOPTIMOD.txt"
    P_SALI = get_Profiles(filename)

    filename=INPUTDIR + "GB_BBP700_MEDSEA_MAY2020_BIOOPTIMOD.txt"
    P_BBP700 = get_Profiles(filename)
    UNIQUE_PROFILES=P_TEMP[:]

  
INDEX_FILE=np.loadtxt(FloatIndexer,dtype=FloatIndex_type, delimiter=",",ndmin=1)



for ip, p in enumerate(UNIQUE_PROFILES):
    for ifile, filename in enumerate(INDEX_FILE['file_name']):
        if (filename.find(p.name()) > -1) & (p.time.strftime("%Y%m%d")==INDEX_FILE[ifile]['time'][:8] ) :
            time= INDEX_FILE[ifile]['time'].replace("-","").replace(":","")
            lon = INDEX_FILE[ifile]['lon'].astype(np.float64)
            lat = INDEX_FILE[ifile]['lat'].astype(np.float64)
            break
    else:
        TI = TimeInterval("20160701","20160811","%Y%m%d")
        if (p.name()=="6901766") & (TI.contains(p.time)):
            print p.ID(), " Ma da dove l'hai preso?"
            continue
        else:
            raise ValueError ( "file %s not found" %(p.name()) )

    wmodir=OUTDIR + os.path.dirname(filename)
    outfile=OUTDIR + filename
    os.system("mkdir -p " + wmodir)
    

    VARLIST=[]

    if setup_2019:
        if p in P_chl  : VARLIST.append('CHL')
        if p in P_Ed380: VARLIST.append('IRR_380')
        if p in P_Ed412: VARLIST.append('IRR_412')
        if p in P_Ed490: VARLIST.append('IRR_490')
        if p in P_PAR: VARLIST.append('PAR')
        same_profile_objects = Same_Profilelist(p, VARLIST, P_chl,P_Ed380, P_Ed412, P_Ed490, P_PAR)


    if setup_2020:
        if p in P_TEMP  : VARLIST.append('TEMP')
        if p in P_SALI  : VARLIST.append('SALI')
        if p in P_BBP700: VARLIST.append('BBP700')
        same_profile_objects = Same_Profilelist_2020(p, VARLIST, P_TEMP, P_SALI, P_BBP700)

    
    V=[]
    for ivar, var in enumerate(VARLIST):
        p = same_profile_objects[ivar]
        Pres, Profile, Qc = p.read(var)
        V.append((Pres,Profile,Qc))

# NetCDF generation
    print outfile
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

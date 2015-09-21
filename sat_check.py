import postproc
from postproc.Timelist import *
import glob,os
import numpy as np
import scipy.io.netcdf as NC


def satread(filename):
    '''
    returns CHL,Lat,Lon
    ''' 
    ncIN = NC.netcdf_file(filename,'r')
    CHL_IN = ncIN.variables['CHL'].data[0,:,:].copy()
    Lon    = ncIN.variables['lon'].data
    Lat    = ncIN.variables['lat'].data
    ncIN.close()
    return CHL_IN,Lat,Lon

def satwrite(filename, CHL,Lon,Lat):
    ncOUT  = NC.netcdf_file(filename,'w')
    ncOUT.createDimension('time', 1)
    ncOUT.createDimension('depth',1)
    ncOUT.createDimension('lon',len(Lon))
    ncOUT.createDimension('lat',len(Lat))
    
    ncvar = ncOUT.createVariable('CHL', 'f', ('lat','lon'))
    ncvar[:] = CHL
    setattr(ncvar, 'missing_value', fillValue)
    setattr(ncvar, '_FillValue', fillValue)
    ncvar = ncOUT.createVariable('depth','f',('depth',))
    ncvar[:] = 1.47210180759
    ncvar = ncOUT.createVariable('time','f',('time',))
    ncvar[:] = 1.0
    ncvar = ncOUT.createVariable('lat','f',('lat',))
    ncvar[:] = Lat
    ncvar = ncOUT.createVariable('lon','f',('lon',))
    ncvar[:] = Lon
    ncOUT.close()        

ORIGDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/ORIG/"
CHECKDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/"
CLIM_FILE="/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/DA/CHECKSAT/SatClimatology.nc"


TL_orig = TimeList("19500101","20500101", ORIGDIR ,"*.nc",'postproc/IOnames_sat.xml')
TLCheck = TimeList("19500101","20500101", CHECKDIR,"*.nc",'postproc/IOnames_sat.xml')

#os.chdir(ORIGDIR); ORIG__LIST=glob.glob("*nc")os.chdir(CHECKDIR); CHECKLIST=glob.glob("*nc")

toCheckList=[]

ORIG_NAMES=[os.path.basename(i) for i in TL_orig.filelist]
CHECKNAMES=[os.path.basename(i) for i in TLCheck.filelist]

for iTimeorig, filename in enumerate(ORIG_NAMES):
    if filename in CHECKNAMES:
        pass
    else:
        toCheckList.append((filename,iTimeorig))

fillValue=-999.0
        
if len(toCheckList)>0:
    ncIN = NC.netcdf_file(CLIM_FILE,'r')
    jpi = ncIN.dimensions['lon']
    jpj = ncIN.dimensions['lat']    
    MEAN = ncIN.variables['Mean'].data #(366, 253, 733)
    STD  = ncIN.variables['Mean'].data 
    ncIN.close()
    
for filename, iTimeorig in toCheckList:
    julian = int( TL_orig.Timelist[iTimeorig].strftime("%j") )
    CHL_IN, Lat,Lon = satread(TL_orig.filelist[iTimeorig])
    CHL_IN[581:,164:] = fillValue # BLACK SEA
    
    DAILY_REF_MEAN = MEAN[julian-1,:,:]
    DAILY_REF_STD  =  STD[julian-1,:,:]
    
    cloudsLandTIME = CHL_IN         == fillValue
    cloudlandsCLIM = DAILY_REF_MEAN == fillValue
    CHL_OUT = CHL_IN.copy()
    CHL_OUT[cloudsLandTIME] = fillValue
    CHL_OUT[cloudlandsCLIM] = fillValue
    counter_refNAN = (~cloudsLandTIME & cloudlandsCLIM).sum(axis=None)
    
    
    outOfRange = np.abs(CHL_IN - DAILY_REF_MEAN) > DAILY_REF_STD *2.0
    outOfRange[cloudsLandTIME | cloudlandsCLIM ] = False
    counter_elim = outOfRange.sum(axis = None)
    CHL_OUT[outOfRange] = fillValue 
    
    print filename
    print 'Rejection:  after check', counter_elim, ' values'
    print 'rejected for NAN in Climatology', counter_refNAN, ' values'
    satwrite(CHECKDIR + filename, CHL_OUT, Lon, Lat)
     
    






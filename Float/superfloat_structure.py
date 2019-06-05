from instruments import bio_float
from instruments import lovbio_float
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import superfloat_generator
import os,sys
import scipy.io.netcdf as NC
import numpy as np
TI = TimeInterval('2012','2020','%Y')
R = Rectangle(-6,36,30,46)

PROFILES_LOV =lovbio_float.FloatSelector('CHLA', TI, R)
OUTDIR="/gpfs/scratch/userexternal/gbolzon0/SuperFloat/" #os.getenv("ONLINE_REPO")


def get_info(p,outdir):
    wmo=p._my_float.wmo
    filename="%s%s/MR%s_%03d.nc" %(outdir,wmo, wmo,p._my_float.cycle)
    return filename

def dumpfile(outfile,p,Pres,chl_profile,Qc):
    PresT, Temp, QcT = p.read('TEMP', read_adjusted=False)
    PresT, Sali, QcS = p.read('PSAL', read_adjusted=False)

    print "dumping " + outfile
    ncOUT = NC.netcdf_file(outfile,"w")
    ncOUT.createDimension("DATETIME",14)
    ncOUT.createDimension("NPROF", 1)
    ncOUT.createDimension('nTEMP', len(PresT))
    ncOUT.createDimension('nPSAL', len(PresT))
    ncOUT.createDimension('nCHLA', len(Pres ))    
    
    ncvar=ncOUT.createVariable("REFERENCE_DATE_TIME", 'c', ("DATETIME",))
    ncvar[:]=p.time.strftime("%Y%m%d%H%M%S")
    ncvar=ncOUT.createVariable("JULD", 'd', ("NPROF",))
    ncvar[:]=0.0
    ncvar=ncOUT.createVariable("LONGITUDE", "d", ("NPROF",))
    ncvar[:] = p.lon.astype(np.float64)
    ncvar=ncOUT.createVariable("LATITUDE", "d", ("NPROF",))
    ncvar[:] = p.lat.astype(np.float64)
 
    
    ncvar=ncOUT.createVariable('TEMP','f',('nTEMP',))
    ncvar[:]=Temp
    ncvar=ncOUT.createVariable('PRES_TEMP','f',('nTEMP',))
    ncvar[:]=PresT
    ncvar=ncOUT.createVariable('TEMP_QC','f',('nTEMP',))
    ncvar[:]=QcT
    
    ncvar=ncOUT.createVariable('PSAL','f',('nTEMP',))
    ncvar[:]=Sali
    ncvar=ncOUT.createVariable('PRES_PSAL','f',('nTEMP',))
    ncvar[:]=PresT
    ncvar=ncOUT.createVariable('PSAL_QC','f',('nTEMP',))
    ncvar[:]=QcS
    
    ncvar=ncOUT.createVariable('CHLA','f',('nCHLA',))
    ncvar[:]=Pres
    ncvar=ncOUT.createVariable('PRES_CHLA','f',('nCHLA',))
    ncvar[:]=chl_profile
    ncvar=ncOUT.createVariable('CHLA_QC','f',('nCHLA',))
    ncvar[:]=Qc
    ncOUT.close()


for ip, pLov in enumerate(PROFILES_LOV[:1]):
    pCor = bio_float.from_lov_profile(pLov, verbose=True)
    is_only_lov = pCor is None
    if is_only_lov:
        outfile = get_info(pLov,OUTDIR)
    else:
        outfile = get_info(pCor,OUTDIR)
    if os.path.exists(outfile): continue
    os.system('mkdir -p ' + os.path.dirname(outfile))
    
    if is_only_lov:
        Pres, Profile, Qc = pLov.read('CHLA',read_adjusted=True)
        profile_for_dump = pLov
    else:
        Pres, Profile,Qc = superfloat_generator.synthesis_profile(pLov, pCor)
        profile_for_dump = pCor
    
    if Pres is None: continue # no data
    
    dumpfile(outfile, profile_for_dump, Pres, Profile, Qc)



# Deve gestire i profili only_coriolis
PROFILES_COR =bio_float.FloatSelector('CHLA', TI, R)

for ip, pCor in enumerate(PROFILES_COR):
    outfile = get_info(pCor, OUTDIR)
    if os.path.exists(outfile): continue
    os.system('mkdir -p ' + os.path.dirname(outfile))
    Pres, CHL, Qc= superfloat_generator.treating_coriolis(pCor)
    if Pres is None: continue # no data
    dumpfile(outfile, pCor, Pres, CHL, Qc)
    
    
    
    

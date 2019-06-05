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
#da rimuovere /gpfs/scratch/userexternal/gbolzon0/SuperFloat/7900592/MR7900592_071.nc

PROFILES_LOV =lovbio_float.FloatSelector('SR_NO3', TI, R)
OUTDIR="/gpfs/scratch/userexternal/gbolzon0/SuperFloat/" #os.getenv("ONLINE_REPO")
def get_info(p,outdir):
    wmo=p._my_float.wmo
    filename="%s%s/MR%s_%03d.nc" %(outdir,wmo, wmo,p._my_float.cycle)
    return filename


def dump_nitrate_file(outfile,pLov, mode='w'):
    Pres, Value, Qc= pLov.read("SR_NO3", read_adjusted=True)
    nP=len(Pres)
    if nP<5 :
        print "few values for " + pLov._my_float.filename
        return
    
    ncOUT = NC.netcdf_file(outfile,mode)
    ncOUT.createDimension('nNITRATE', nP)
    ncvar=ncOUT.createVariable("PRES_NITRATE", 'f', ('nNITRATE',))
    ncvar[:]=Pres
    ncvar=ncOUT.createVariable("NITRATE", 'f', ('nNITRATE',))
    ncvar[:]=Value
    ncvar=ncOUT.createVariable("NITRATE", 'f', ('nNITRATE',))
    ncvar[:]=Qc
    ncOUT.close()
def exist_nitrate(filename):
    ncIN = NC.netcdf_file(filename,'r')
    variables=ncIN.variables.keys()
    ncIN.close()
    return 'NITRATE' in variables

force_writing_nitrate=False

for ip, pLov in enumerate(PROFILES_LOV):
    pCor = bio_float.from_lov_profile(pLov, verbose=False)
    is_only_lov = pCor is None
    if is_only_lov: # wmo = [6903235, 6902900 6902902]
        outfile = get_info(pLov,OUTDIR)
    else:
        outfile = get_info(pCor,OUTDIR)
    if os.path.exists(outfile):
        if not exist_nitrate(outfile):
            dumpfile(outfile, pLov, mode='a')
        else:
            if force_writing_nitrate:
                dump_nitrate_file(outfile, pLov, mode='a')
    else:
        print outfile + " not found"
        dump_nitrate_file(outfile, pLov, mode='w')
    
        
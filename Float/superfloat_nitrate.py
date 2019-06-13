from instruments import bio_float
from instruments import lovbio_float
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import superfloat_generator
import os,sys
import scipy.io.netcdf as NC
import numpy as np

def get_outfile(p,outdir):
    wmo=p._my_float.wmo
    filename="%s%s/MR%s_%03d.nc" %(outdir,wmo, wmo,p._my_float.cycle)
    return filename


def dump_nitrate_file(outfile, p_pos, p, Pres, Value, Qc, metadata,mode='w'):
    
    nP=len(Pres)
    ncOUT = NC.netcdf_file(outfile,mode)
    if mode=='w': # if not existing file, we'll put header, TEMP, PSAL
        PresT, Temp, QcT = p.read('TEMP', read_adjusted=False)
        PresT, Sali, QcS = p.read('PSAL', read_adjusted=False)
        ncOUT.createDimension("DATETIME",14)
        ncOUT.createDimension("NPROF", 1)
        ncOUT.createDimension('nTEMP', len(PresT))
        ncOUT.createDimension('nPSAL', len(PresT))

        ncvar=ncOUT.createVariable("REFERENCE_DATE_TIME", 'c', ("DATETIME",))
        ncvar[:]=p_pos.time.strftime("%Y%m%d%H%M%S")
        ncvar=ncOUT.createVariable("JULD", 'd', ("NPROF",))
        ncvar[:]=0.0
        ncvar=ncOUT.createVariable("LONGITUDE", "d", ("NPROF",))
        ncvar[:] = p_pos.lon.astype(np.float64)
        ncvar=ncOUT.createVariable("LATITUDE", "d", ("NPROF",))
        ncvar[:] = p_pos.lat.astype(np.float64)


        ncvar=ncOUT.createVariable('TEMP','f',('nTEMP',))
        ncvar[:]=Temp
        setattr(ncvar, 'origin'     , metadata.origin)
        setattr(ncvar, 'file_origin', metadata.filename)
        setattr(ncvar, 'variable'   , 'TEMP')
        setattr(ncvar, 'units'      , "degree_Celsius")
        ncvar=ncOUT.createVariable('PRES_TEMP','f',('nTEMP',))
        ncvar[:]=PresT
        ncvar=ncOUT.createVariable('TEMP_QC','f',('nTEMP',))
        ncvar[:]=QcT

        ncvar=ncOUT.createVariable('PSAL','f',('nTEMP',))
        ncvar[:]=Sali
        setattr(ncvar, 'origin'     , metadata.origin)
        setattr(ncvar, 'file_origin', metadata.filename)
        setattr(ncvar, 'variable'   , 'SALI')
        setattr(ncvar, 'units'      , "PSS78")
        ncvar=ncOUT.createVariable('PRES_PSAL','f',('nTEMP',))
        ncvar[:]=PresT
        ncvar=ncOUT.createVariable('PSAL_QC','f',('nTEMP',))
        ncvar[:]=QcS


    ncOUT.createDimension('nNITRATE', nP)
    ncvar=ncOUT.createVariable("PRES_NITRATE", 'f', ('nNITRATE',))
    ncvar[:]=Pres
    ncvar=ncOUT.createVariable("NITRATE", 'f', ('nNITRATE',))
    ncvar[:]=Value
    setattr(ncvar, 'origin'     , metadata.origin)
    setattr(ncvar, 'file_origin', metadata.filename)
    setattr(ncvar, 'variable'   , 'SR_NO3_ADJUSTED')
    setattr(ncvar, 'units'      , "mmol/m3")
    ncvar=ncOUT.createVariable("NITRATE_QC", 'f', ('nNITRATE',))
    ncvar[:]=Qc
    ncOUT.close()


TI = TimeInterval('2012','2020','%Y')
R = Rectangle(-6,36,30,46)


PROFILES_LOV =lovbio_float.FloatSelector('SR_NO3', TI, R)
wmo_list= lovbio_float.get_wmo_list(PROFILES_LOV)

OUTDIR="/gpfs/scratch/userexternal/gbolzon0/SuperFloat/" #os.getenv("ONLINE_REPO")

force_writing_nitrate=True

for wmo in wmo_list:
    Profilelist = lovbio_float.filter_by_wmo(PROFILES_LOV, wmo)
    for ip, pLov in enumerate(Profilelist):
        pCor = bio_float.from_lov_profile(pLov, verbose=True)
        is_only_lov = pCor is None
        if is_only_lov: # wmo = [6903235, 6902900 6902902]
            outfile = get_outfile(pLov,OUTDIR)
            profile_for_position=pLov
            #print "from LOV ", outfile
        else:
            outfile = get_outfile(pCor,OUTDIR)
            profile_for_position=pCor
            #print "from COR ", outfile
        metatata = superfloat_generator.Metadata('lov',pLov._my_float.filename)

        Pres, Value, Qc= pLov.read("SR_NO3", read_adjusted=True)
        nP=len(Pres)
        if nP<5 :
            print "few values for " + pLov._my_float.filename
            continue
        os.system('mkdir -p ' + os.path.dirname(outfile))

        if os.path.exists(outfile):
            if not superfloat_generator.exist_variable('NITRATE', outfile):
                dump_nitrate_file(outfile, profile_for_position, pLov, Pres, Value, Qc, metatata,mode='a')
            else:
                if force_writing_nitrate:
                    dump_nitrate_file(outfile, profile_for_position, pLov, Pres, Value, Qc, metatata,mode='a')
        else:
            print outfile + " not found"
            dump_nitrate_file(outfile, profile_for_position, pLov,Pres, Value, Qc, metatata,mode='w')


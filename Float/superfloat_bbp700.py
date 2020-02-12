import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates superfloat files of BBP700.
    Reads from Coriolis.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--datestart','-s',
                                type = str,
                                required = True,
                                help = '''date in yyyymmdd format''')
    parser.add_argument(   '--dateend','-e',
                                type = str,
                                required = True,
                                help = '''date in yyyymmdd format ''')
    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                default = "/gpfs/scratch/userexternal/gbolzon0/SUPERFLOAT/",
                                help = 'path of the Superfloat dataset ')

    return parser.parse_args()

args = argument()

from instruments import bio_float
from instruments import lovbio_float
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import superfloat_generator
from commons.utils import addsep
import os
import scipy.io.netcdf as NC
import numpy as np
import seawater as sw


def dump_bbp700_file(outfile, p, Pres, Value, Qc, metatata, mode='w'):
    nP=len(Pres)
    if mode=='a':
        command = "cp %s %s.tmp" %(outfile,outfile)
        os.system(command)
    ncOUT = NC.netcdf_file(outfile + ".tmp" ,mode)

    if mode=='w': # if not existing file, we'll put header, TEMP, PSAL
        PresT, Temp, QcT = p.read('TEMP', read_adjusted=False)
        PresT, Sali, QcS = p.read('PSAL', read_adjusted=False)        
        ncOUT.createDimension("DATETIME",14)
        ncOUT.createDimension("NPROF", 1)
        ncOUT.createDimension('nTEMP', len(PresT))
        ncOUT.createDimension('nPSAL', len(PresT))

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

    print "dumping bbp700 on " + outfile
    bbp700_already_existing="nBBP700" in ncOUT.dimensions.keys()
    if not bbp700_already_existing : ncOUT.createDimension('nBBP700', nP)
    ncvar=ncOUT.createVariable("PRES_BBP700", 'f', ('nBBP700',))
    ncvar[:]=Pres
    ncvar=ncOUT.createVariable("BBP700", 'f', ('nBBP700',))
    ncvar[:]=Value
    if not bbp700_already_existing:
        setattr(ncvar, 'origin'     , metadata.origin)
        setattr(ncvar, 'file_origin', metadata.filename)
        setattr(ncvar, 'variable'   , 'BBP700')
        setattr(ncvar, 'units'      , "m-1")
        setattr(ncvar, 'longname'   , 'Particle backscattering at 700 nanometers')
    ncvar=ncOUT.createVariable("BBP700_QC", 'f', ('nBBP700',))
    ncvar[:]=Qc
    ncOUT.close()

    os.system("mv " + outfile + ".tmp " + outfile)


OUTDIR = addsep(args.outdir)
TI     = TimeInterval(args.datestart,args.dateend,'%Y%m%d')
R = Rectangle(-6,36,30,46)
force_writing_=False

PROFILES_COR =bio_float.FloatSelector('BBP700', TI, R)

wmo_list= lovbio_float.get_wmo_list(PROFILES_COR)



def get_outfile(p,outdir):
    wmo=p._my_float.wmo
    filename="%s%s/MR%s_%03d.nc" %(outdir,wmo, wmo,p._my_float.cycle)
    return filename

for wmo in wmo_list:
    print wmo
    Profilelist=bio_float.filter_by_wmo(PROFILES_COR, wmo)
    for ip, pCor in enumerate(Profilelist):
        outfile = get_outfile(pCor,OUTDIR)
        metadata = superfloat_generator.Metadata('Coriolis', pCor._my_float.filename)
        os.system('mkdir -p ' + os.path.dirname(outfile))

        if superfloat_generator.exist_valid(outfile):
            if not superfloat_generator.exist_variable('BBP700', outfile):
                Pres, Value, Qc = pCor.read('BBP700', read_adjusted=False)
                if Pres is not None: dump_bbp700_file(outfile, pCor, Pres, Value, Qc, metadata,mode='a')
            else:
                if force_writing_:
                    Pres, Value, Qc = pCor.read('BBP700', read_adjusted=False)
                    if Pres is not None: dump_bbp700_file(outfile, pCor, Pres, Value, Qc, metadata,mode='a')
        else:
            Pres, Value, Qc = pCor.read('BBP700', read_adjusted=False)
            if Pres is not None: dump_bbp700_file(outfile, pCor, Pres, Value, Qc, metadata,mode='w')


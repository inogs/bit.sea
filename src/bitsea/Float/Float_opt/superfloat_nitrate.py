import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates superfloat files of chla.
    Reads from LOV dataset.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--datestart','-s',
                                type = str,
                                required = True,
                                help = '''date in yyyymmdd format''')
    parser.add_argument(   '--dateend','-e',
                                type = str,
                                required = True,
                                help = '''date in yyyymmdd format''')
    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                default = "/gpfs/scratch/userexternal/gbolzon0/SUPERFLOAT/",
                                help = 'path of the Superfloat dataset ')
    parser.add_argument(   '--force', '-f',
                                action='store_true',
                                help = """Overwrite existing files
                                """)

    return parser.parse_args()

args = argument()

from instruments import optbio_float_2020
from instruments import lovbio_float
from commons.time_interval import TimeInterval
from basins.region import Rectangle
from Float import superfloat_generator
from commons.utils import addsep
import os
import scipy.io.netcdf as NC
import numpy as np

def get_outfile(p,outdir):
    filename="%s%s" %(outdir,p._my_float.filename)
    return filename


def dump_nitrate_file(outfile, p_pos, p, Pres, Value, Qc, metadata,mode='w'):
    
    nP=len(Pres)
    if mode=='a':
        command = "cp %s %s.tmp" %(outfile,outfile)
        os.system(command)
    ncOUT = NC.netcdf_file(outfile + ".tmp",mode)
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

    print "dumping nitrate on " + outfile
    nitrate_already_existing="nNITRATE" in ncOUT.dimensions.keys()
    if not nitrate_already_existing: ncOUT.createDimension('nNITRATE', nP)
    ncvar=ncOUT.createVariable("PRES_NITRATE", 'f', ('nNITRATE',))
    ncvar[:]=Pres
    ncvar=ncOUT.createVariable("NITRATE", 'f', ('nNITRATE',))
    ncvar[:]=Value
    if not nitrate_already_existing:
        setattr(ncvar, 'origin'     , metadata.origin)
        setattr(ncvar, 'file_origin', metadata.filename)
        setattr(ncvar, 'variable'   , 'SR_NO3_ADJUSTED')
        setattr(ncvar, 'units'      , "mmol/m3")
    ncvar=ncOUT.createVariable("NITRATE_QC", 'f', ('nNITRATE',))
    ncvar[:]=Qc
    ncOUT.close()
    os.system("mv " + outfile + ".tmp " + outfile)

OUTDIR = addsep(args.outdir)
TI     = TimeInterval(args.datestart,args.dateend,'%Y%m%d')
R = Rectangle(-6,36,30,46)


PROFILES_LOV =lovbio_float.FloatSelector('SR_NO3', TI, R)
wmo_list= lovbio_float.get_wmo_list(PROFILES_LOV)
nWMOS=len(wmo_list)

force_writing_nitrate=args.force

for iwmo, wmo in enumerate(wmo_list):
    print wmo, iwmo, " of ", nWMOS
    Profilelist = lovbio_float.filter_by_wmo(PROFILES_LOV, wmo)
    for ip, pLov in enumerate(Profilelist):
        p_2020 = optbio_float_2020.from_profile(pLov, verbose=True)
        if p_2020 is None: continue
        is_only_lov = p_2020 is None
        outfile = get_outfile(p_2020,OUTDIR)

        metatata = superfloat_generator.Metadata('lov',pLov._my_float.filename)

        Pres, Value, Qc= pLov.read("SR_NO3", read_adjusted=True)
        nP=len(Pres)
        if nP<5 :
            print "few values for " + pLov._my_float.filename
            continue
        #os.system('mkdir -p ' + os.path.dirname(outfile))

        if superfloat_generator.exist_valid(outfile):
            if not superfloat_generator.exist_variable('NITRATE', outfile):
                dump_nitrate_file(outfile, p_2020, pLov, Pres, Value, Qc, metatata,mode='a')
            else:
                if force_writing_nitrate:
                    dump_nitrate_file(outfile, p_2020, pLov, Pres, Value, Qc, metatata,mode='a')
#         else:
#             print outfile + " not found"
#             dump_nitrate_file(outfile, profile_for_position, pLov,Pres, Value, Qc, metatata,mode='w')


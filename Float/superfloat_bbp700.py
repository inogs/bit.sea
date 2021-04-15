import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates superfloat files of BBP700.
    Reads from Coriolis.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--datestart','-s',
                                type = str,
                                required = False,
                                help = '''date in yyyymmdd format''')
    parser.add_argument(   '--dateend','-e',
                                type = str,
                                required = False,
                                help = '''date in yyyymmdd format ''')
    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                default = "/gpfs/scratch/userexternal/gbolzon0/SUPERFLOAT/",
                                help = 'path of the Superfloat dataset ')
    parser.add_argument(   '--force', '-f',
                                action='store_true',
                                help = """Overwrite existing files
                                """)
    parser.add_argument(   '--update_file','-u',
                                type = str,
                                required = False,
                                default = 'NO_file',
                                help = '''file with updated floats''')
    return parser.parse_args()

args = argument()

if (args.datestart == 'NO_data') & (args.dateend == 'NO_data') & (args.update_file == 'NO_file'):
    raise ValueError("No file nor data inserted: you have to pass either datastart and dataeend or the update_file")

if ((args.datestart == 'NO_data') or (args.dateend == 'NO_data')) & (args.update_file == 'NO_file'):
    raise ValueError("No file nor data inserted: you have to pass both datastart and dataeend")

from instruments import bio_float
#from instruments import lovbio_float
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import superfloat_generator
from commons.utils import addsep
import os
import scipy.io.netcdf as NC
import numpy as np
import datetime


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

def get_outfile(p,outdir):
    wmo=p._my_float.wmo
    filename="%s%s/%s" %(outdir,wmo, os.path.basename(p._my_float.filename))
    return filename

input_file=args.update_file
if input_file == 'NO_file':

    OUTDIR = addsep(args.outdir)
    TI     = TimeInterval(args.datestart,args.dateend,'%Y%m%d')
    R = Rectangle(-6,36,30,46)
    force_writing_bbp=args.force

    PROFILES_COR =bio_float.FloatSelector('BBP700', TI, R)

    wmo_list= bio_float.get_wmo_list(PROFILES_COR)



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
                    if force_writing_bbp:
                        Pres, Value, Qc = pCor.read('BBP700', read_adjusted=False)
                        if Pres is not None: dump_bbp700_file(outfile, pCor, Pres, Value, Qc, metadata,mode='a')
            else:
                Pres, Value, Qc = pCor.read('BBP700', read_adjusted=False)
                if Pres is not None: dump_bbp700_file(outfile, pCor, Pres, Value, Qc, metadata,mode='w')

else:

    OUTDIR = addsep(args.outdir)
    mydtype= np.dtype([
        ('file_name','S200'),
        ('date','S200'),
        ('latitude',np.float32),
        ('longitude',np.float32),
        ('ocean','S10'),
        ('profiler_type',np.int),
        ('institution','S10'),
        ('parameters','S200'),
        ('parameter_data_mode','S100'),
        ('date_update','S200')] )

    INDEX_FILE=np.loadtxt(input_file,dtype=mydtype, delimiter=",",ndmin=1,skiprows=0)
    nFiles=INDEX_FILE.size

    for iFile in range(nFiles):
        timestr          = INDEX_FILE['date'][iFile]
        lon              = INDEX_FILE['longitude' ][iFile]
        lat              = INDEX_FILE['latitude' ][iFile]
        filename         = INDEX_FILE['file_name'][iFile]
        available_params = INDEX_FILE['parameters'][iFile]
        parameterdatamode= INDEX_FILE['parameter_data_mode'][iFile]
        float_time = datetime.datetime.strptime(timestr,'%Y%m%d%H%M%S')
        filename=filename.replace('coriolis/','').replace('profiles/','')

        if  'BBP700' in available_params:
            pCor=bio_float.profile_gen(lon, lat, float_time, filename, available_params,parameterdatamode)
            outfile = get_outfile(pCor,OUTDIR)
            os.system('mkdir -p ' + os.path.dirname(outfile))
            if pCor._my_float.status_var('BBP700') in ['A', 'D']:
                Pres, Value, Qc = pCor.read('BBP700', read_adjusted=True)
            else:
                Pres, Value, Qc = pCor.read('BBP700', read_adjusted=False)
            metadata = superfloat_generator.Metadata('Coriolis', pCor._my_float.filename)
            if Pres is None: continue # no data

            writing_mode='w'
            if superfloat_generator.exist_valid(outfile): writing_mode='a'
            dump_bbp700_file(outfile, pCor, Pres, Value, Qc, metadata,mode=writing_mode)


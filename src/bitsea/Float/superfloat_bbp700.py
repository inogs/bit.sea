import argparse
from bitsea.utilities.argparse_types import existing_dir_path
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
                                type = existing_dir_path,
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

from bitsea.instruments import bio_float
from bitsea.commons.time_interval import TimeInterval
from bitsea.basins.region import Rectangle
from bitsea.Float import superfloat_generator
from pathlib import Path
import os
import netCDF4 as NC
import numpy as np
import datetime

class Metadata():
    def __init__(self, filename):
        self.filename = str(filename)
        self.status_var = 'n'

def dump_bbp700_file(outfile, p, Pres, Value, Qc, metadata, mode='w'):
    nP=len(Pres)
    tmpfile=outfile.with_suffix('.nc.tmp')
    if mode=='a':
        command = "cp {} {}".format(outfile,tmpfile)
        os.system(command)
    ncOUT = NC.Dataset(tmpfile ,mode)

    if mode=='w':# if not existing file, we'll put header, TEMP, PSAL
        setattr(ncOUT, 'origin'     , 'coriolis')
        setattr(ncOUT, 'file_origin', metadata.filename)
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
        setattr(ncvar, 'variable'   , 'TEMP')
        setattr(ncvar, 'units'      , "degree_Celsius")
        ncvar=ncOUT.createVariable('PRES_TEMP','f',('nTEMP',))
        ncvar[:]=PresT
        ncvar=ncOUT.createVariable('TEMP_QC','f',('nTEMP',))
        ncvar[:]=QcT

        ncvar=ncOUT.createVariable('PSAL','f',('nTEMP',))
        ncvar[:]=Sali
        setattr(ncvar, 'variable'   , 'SALI')
        setattr(ncvar, 'units'      , "PSS78")
        ncvar=ncOUT.createVariable('PRES_PSAL','f',('nTEMP',))
        ncvar[:]=PresT
        ncvar=ncOUT.createVariable('PSAL_QC','f',('nTEMP',))
        ncvar[:]=QcS

    print("dumping bbp700 on " + str(outfile), flush=True)
    bbp700_already_existing=superfloat_generator.exist_valid_variable("BBP700", tmpfile)
    if not bbp700_already_existing :
        ncOUT.createDimension('nBBP700', nP)
        ncvar=ncOUT.createVariable("PRES_BBP700", 'f', ('nBBP700',))
        ncvar[:]=Pres
        ncvar=ncOUT.createVariable("BBP700", 'f', ('nBBP700',))
        ncvar[:]=Value
        setattr(ncvar, 'status_var' , metadata.status_var)
        setattr(ncvar, 'variable'   , 'BBP700')
        setattr(ncvar, 'units'      , "m-1")
        setattr(ncvar, 'longname'   , 'Particle backscattering at 700 nanometers')
        ncvar=ncOUT.createVariable("BBP700_QC", 'f', ('nBBP700',))
        ncvar[:]=Qc
    else:
        ncvar=ncOUT.variables['PRES_BBP700']
        ncvar[:]=Pres
        ncvar=ncOUT.variables['BBP700']
        ncvar[:]=Value
        ncvar=ncOUT.variables['BBP700_QC']
        ncvar[:]=Qc

    ncOUT.close()

    os.system("mv {} {}".format(tmpfile,outfile))

def get_outfile(p,outdir):
    wmo=p._my_float.wmo
    filename = outdir / wmo / p._my_float.filename.name
    return filename

def bbp_algorithm(pCor, outfile, metadata,writing_mode):
    Pres, _, _ = pCor.read('TEMP', read_adjusted=False)
    if len(Pres)<5:
        print("few values in Coriolis TEMP in " + str(pCor._my_float.filename), flush=True)
        return
    Path.mkdir(outfile.parent, exist_ok=True)
    metadata.status_var = pCor._my_float.status_var('BBP700')
    if metadata.status_var in ['A', 'D']:
        Pres, Value, Qc = pCor.read('BBP700', read_adjusted=True)
    else:
        Pres, Value, Qc = pCor.read('BBP700', read_adjusted=False)
    if Pres is None: return
    if len(Pres)<5:
        print("few values in Coriolis for BBP in " + str(pCor._my_float.filename), flush=True)
        return
    dump_bbp700_file(outfile, pCor, Pres, Value, Qc, metadata,mode=writing_mode)

OUTDIR = args.outdir
input_file=args.update_file

if input_file == 'NO_file':

    TI     = TimeInterval(args.datestart,args.dateend,'%Y%m%d')
    R = Rectangle(-6,36,30,46)

    PROFILES_COR =bio_float.FloatSelector('BBP700', TI, R)

    wmo_list= bio_float.get_wmo_list(PROFILES_COR)
    wmo_list.sort()

    for wmo in wmo_list:
        print(wmo, flush=True)
        Profilelist=bio_float.filter_by_wmo(PROFILES_COR, wmo)
        for ip, pCor in enumerate(Profilelist):
            outfile = get_outfile(pCor,OUTDIR)
            writing_mode=superfloat_generator.writing_mode(outfile)

            condition_to_write = ~superfloat_generator.exist_valid_variable('BBP700',outfile)
            if args.force: condition_to_write=True
            if not condition_to_write: continue

            metadata = Metadata(pCor._my_float.filename)
            bbp_algorithm(pCor, outfile, metadata, writing_mode)

else:

    INDEX_FILE=superfloat_generator.read_float_update(input_file)

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

        if 'BBP700' in available_params:
            pCor=bio_float.profile_gen(lon, lat, float_time, filename, available_params,parameterdatamode)
            outfile = get_outfile(pCor,OUTDIR)
            writing_mode=superfloat_generator.writing_mode(outfile)

            metadata = Metadata(pCor._my_float.filename)
            bbp_algorithm(pCor, outfile, metadata, writing_mode)

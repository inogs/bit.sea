from bitsea.utilities.argparse_types import existing_dir_path
import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates superfloat files of chla.
    Reads from Coriolis dataset.
    It has to be called as the first one of the series of superfloat generators.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--datestart','-s',
                                type = str,
                                required = False,
                                default = 'NO_data',
                                help = '''date in yyyymmss format''')
    parser.add_argument(   '--dateend','-e',
                                type = str,
                                required = False,
                                default = 'NO_data',
                                help = '''date in yyyymmss format''')
    parser.add_argument(   '--outdir','-o',
                                type = existing_dir_path,
                                required = True,                                
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


from pathlib import Path
from bitsea.instruments import bio_float
from bitsea.commons.time_interval import TimeInterval
from bitsea.basins.region import Rectangle
from bitsea.Float import superfloat_generator
import netCDF4 as NC
import numpy as np
import datetime

class Metadata():
    def __init__(self, filename):
        self.filename = str(filename)
        self.status_var = 'n'
        self.shift_from_coriolis=0

def get_outfile(p,outdir):
    wmo=p._my_float.wmo
    filename = outdir / wmo / p._my_float.filename.name
    return filename

def dumpfile(outfile, p,Pres,chl_profile,Qc,metadata):
    _    , Temp, QcT = p.read('TEMP', read_adjusted=False)
    PresT, Sali, QcS = p.read('PSAL', read_adjusted=False)

    print("dumping chla on " , outfile, p.time.strftime(" %Y%m%d-%H:%M:%S"), flush=True)
    ncOUT = NC.Dataset(outfile,"w")
    setattr(ncOUT, 'origin'     , 'coriolis')
    setattr(ncOUT, 'file_origin', metadata.filename)

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
    ncvar[:] = Temp
    setattr(ncvar, 'variable'   , 'TEMP')
    setattr(ncvar, 'units'      , "degree_Celsius")

    ncvar=ncOUT.createVariable('PRES_TEMP','f',('nTEMP',))
    ncvar[:]=PresT
    ncvar=ncOUT.createVariable('TEMP_QC','f',('nTEMP',))
    ncvar[:]=QcT
    
    ncvar=ncOUT.createVariable('PSAL','f',('nPSAL',))
    ncvar[:]=Sali
    setattr(ncvar, 'variable'   , 'SALI')
    setattr(ncvar, 'units'      , "PSS78")

    ncvar=ncOUT.createVariable('PRES_PSAL','f',('nPSAL',))
    ncvar[:]=PresT
    ncvar=ncOUT.createVariable('PSAL_QC','f',('nPSAL',))
    ncvar[:]=QcS

    ncvar=ncOUT.createVariable('CHLA','f',('nCHLA',))
    ncvar[:]=chl_profile
    setattr(ncvar, 'status_var' , metadata.status_var)
    setattr(ncvar, 'shift_from_coriolis', str(metadata.shift_from_coriolis))
    setattr(ncvar, 'variable'   , 'CHLA_ADJUSTED')
    setattr(ncvar, 'units'      , "milligram/m3")

    ncvar=ncOUT.createVariable('PRES_CHLA','f',('nCHLA',))
    ncvar[:]=Pres
    ncvar=ncOUT.createVariable('CHLA_QC','f',('nCHLA',))
    ncvar[:]=Qc
    ncOUT.close()


def chla_algorithm(pCor,outfile):
    '''
    Rejecting profiles with
    - less than 5 CHLA values
    - less than 5 TEMP values
    - min of depth > 200m
    - RealTime Mode

    Shifts the profile in order to have mean=0 in 400-600m
    Negative values are replaced by a background value of 0.005 mg/m3
    '''

    background_value=0.005
    Pres, _, _ = pCor.read('TEMP', read_adjusted=False)
    if len(Pres)<5:
        print("few values in Coriolis TEMP in " + str(pCor._my_float.filename), flush=True)
        return

    metadata = Metadata(pCor._my_float.filename)
    metadata.status_var = pCor._my_float.status_var('CHLA')
    if pCor._my_float.status_var('CHLA') in ['A','D'] :
        Pres,CHL, Qc=pCor.read('CHLA', read_adjusted=True)
        if len(Pres)<5:
            print("few values in Coriolis in CHLA for " + str(pCor._my_float.filename), flush=True)
            return
        if Pres.min() > 200:
            print("Only deep values in Coriolis in CHLA for " + str(pCor._my_float.filename), flush=True)
            return

        ii=(Pres >= 400) & (Pres <= 600)
        if ii.sum() > 0:
            shift = CHL[ii].mean()
            metadata.shift_from_coriolis = shift
            CHL = CHL - shift
        ii=CHL<=0
        CHL[ii] = background_value
    else:
        print("R -- not dumped ", pCor._my_float.filename, flush=True)


    Path.mkdir(outfile.parent, exist_ok=True)
    dumpfile(outfile, pCor, Pres, CHL, Qc, metadata)

OUTDIR = args.outdir
input_file=args.update_file

if input_file == 'NO_file': 
    TI = TimeInterval(args.datestart,args.dateend,'%Y%m%d')
    R = Rectangle(-6,36,30,46)
    PROFILES_COR =bio_float.FloatSelector('CHLA', TI, R)
    wmo_list= bio_float.get_wmo_list(PROFILES_COR)
    wmo_list.sort()

    for wmo in wmo_list:
        Profilelist = bio_float.filter_by_wmo(PROFILES_COR, wmo)
        for ip, pCor in enumerate(Profilelist):
            outfile = get_outfile(pCor, OUTDIR)
            if superfloat_generator.exist_valid(outfile) & (not args.force): continue
            chla_algorithm(pCor,outfile)


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
        
        if 'CHLA' in available_params:
            pCor=bio_float.profile_gen(lon, lat, float_time, filename, available_params,parameterdatamode)
            outfile = get_outfile(pCor, OUTDIR)
            chla_algorithm(pCor, outfile)


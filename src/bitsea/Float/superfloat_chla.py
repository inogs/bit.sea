from bitsea.utilities.argparse_types import existing_dir_path
import argparse
import os
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd


def read_temp_psal(p):
    PresT, Temp, QcT = p.read('TEMP', read_adjusted=True)
    Pres, Sali, QcS = p.read('PSAL', read_adjusted=True)
    if (Pres is None or PresT is None or Temp is None or Sali is None or len(Pres) < 5 or len(PresT) < 5):
        PresT, Temp, QcT = p.read('TEMP', read_adjusted=False)
        Pres, Sali, QcS = p.read('PSAL', read_adjusted=False)
    return PresT, Temp, QcT, Pres, Sali, QcS

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
    parser.add_argument(   '--outdiag','-O',
                                type = existing_dir_path,
                                required = True,
                                help = 'path for statistics, diagnostics, logs')
    parser.add_argument(   '--force', '-f',
                                action='store_true',
                                help = """Overwrite existing files
                                """)

    parser.add_argument(   '--update_file','-u',
                                type = str,
                                required = False,
                                default = 'NO_file',
                                help = '''file with updated floats''')
    parser.add_argument(   '--workers','-w',
                                type = int,
                                required = False,
                                default = os.cpu_count() or 1,
                                help = 'number of worker processes (one float, i.e. WMO, per task)')

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

CHLA_MAX_THRESHOLD = 7.0  # profiles with any CHLA value above this are rejected

class Metadata():
    def __init__(self, filename):
        self.filename = str(filename)
        self.status_var = 'n'
        self.shift_from_coriolis=0

def get_outfile(p,outdir):
    wmo=p._my_float.wmo
    filename = outdir / wmo / p._my_float.filename.name
    return filename


def save_report(out_dir, df_report, filename='Floats_chla_report.csv'):
    """Append the rows of df_report to a CSV, avoiding duplicate rows.
    Called only by the main process, so the CSV has a single writer."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    file_path = out_dir / filename

    df_new = df_report.copy()
    if file_path.exists():
        df_existing = pd.read_csv(file_path)
        df_new = pd.concat([df_existing, df_new], ignore_index=True)

    df_new.drop_duplicates(inplace=True)
    df_new.to_csv(file_path, index=False)


def build_report_row(p, outcome, reason='', shift=np.nan, replaced_values=0):
    return pd.DataFrame([{
        'WMO': p._my_float.wmo,
        'cycle': p._my_float.cycle,
        'time': str(p.time),
        'file_origin': str(p._my_float.filename),
        'outcome': outcome,
        'reason': reason,
        'shift_from_coriolis': shift,
        'background_value': 0.005 if replaced_values else np.nan,
        'replaced_values': replaced_values,
    }])


def dumpfile(outfile, p,Pres,chl_profile,Qc,metadata):
    PresT, Temp, QcT, PresS, Sali, QcS = read_temp_psal(p)

    print("dumping chla on " , outfile, p.time.strftime(" %Y%m%d-%H:%M:%S"), flush=True)
    ncOUT = NC.Dataset(outfile,"w")
    setattr(ncOUT, 'origin'     , 'coriolis')
    setattr(ncOUT, 'file_origin', metadata.filename)

    ncOUT.createDimension("DATETIME",14)
    ncOUT.createDimension("NPROF", 1)

    ncOUT.createDimension('nTEMP', len(PresT))
    ncOUT.createDimension('nPSAL', len(PresS))
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
    ncvar[:]=PresS
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
    - at least one CHLA value > CHLA_MAX_THRESHOLD
    - RealTime Mode

    Shifts the profile in order to have mean=0 in 400-600m
    Negative values are replaced by a background value of 0.005 mg/m3
    '''

    background_value=0.005
    PresT, Temp, _, _, _, _ = read_temp_psal(pCor)
    if PresT is None or Temp is None or len(PresT) < 5 or len(Temp) < 5:
        print("few values in Coriolis TEMP in " + str(pCor._my_float.filename), flush=True)
        return build_report_row(pCor, 'rejected', reason='insufficient_TEMP_data')

    metadata = Metadata(pCor._my_float.filename)
    metadata.status_var = pCor._my_float.status_var('CHLA')
    if metadata.status_var in ['A','D'] :
        Pres,CHL, Qc=pCor.read('CHLA', read_adjusted=True)
        if Pres is None or CHL is None or Qc is None or len(Pres)<5 or len(CHL)<5:
            print("few values in Coriolis in CHLA for " + str(pCor._my_float.filename), flush=True)
            return build_report_row(pCor, 'rejected', reason='insufficient_CHLA_data')
        if len(Pres) != len(CHL) or len(Pres) != len(Qc):
            return build_report_row(pCor, 'rejected', reason='inconsistent_CHLA_array_lengths')
        if Pres.min() > 200:
            print("Only deep values in Coriolis in CHLA for " + str(pCor._my_float.filename), flush=True)
            return build_report_row(pCor, 'rejected', reason='CHLA_only_below_200m')
        if np.any(CHL > CHLA_MAX_THRESHOLD):
            print("CHLA value above %s in %s" % (CHLA_MAX_THRESHOLD, pCor._my_float.filename), flush=True)
            return build_report_row(pCor, 'rejected', reason='CHLA_above_max_threshold')

        ii=(Pres >= 400) & (Pres <= 600)
        shift = 0.0
        if ii.sum() > 0:
            shift = CHL[ii].mean()
            metadata.shift_from_coriolis = shift
            CHL = CHL - shift
        ii=CHL<=0
        replaced_values = int(ii.sum())
        CHL[ii] = background_value
    else:
        print("R -- not dumped ", pCor._my_float.filename, flush=True)
        return build_report_row(pCor, 'rejected', reason=f"unsupported_CHLA_status_{metadata.status_var}")


    Path.mkdir(outfile.parent, exist_ok=True)
    dumpfile(outfile, pCor, Pres, CHL, Qc, metadata)
    return build_report_row(
        pCor, 'accepted', shift=metadata.shift_from_coriolis,
        replaced_values=replaced_values)

OUTDIR = args.outdir
OUTDIAG = args.outdiag
input_file=args.update_file

# wmo -> list of profiles. Filled before the pool is created, so that the
# forked workers inherit it and only the wmo is sent to them.
WMO_PROFILES = {}

def process_wmo(wmo):
    """Processes all the profiles of one float. Returns the report rows;
    the CSV is written only by the main process."""
    rows = []
    for pCor in WMO_PROFILES[wmo]:
        outfile = get_outfile(pCor, OUTDIR)
        if input_file == 'NO_file' and superfloat_generator.exist_valid(outfile) and not args.force:
            continue
        try:
            rows.append(chla_algorithm(pCor, outfile))
        except Exception as exc:
            print("error on " + str(pCor._my_float.filename) + ": " + repr(exc), flush=True)
            rows.append(build_report_row(pCor, 'error', reason=repr(exc)))
    return rows

if input_file == 'NO_file': 
    TI = TimeInterval(args.datestart,args.dateend,'%Y%m%d')
    R = Rectangle(-6,36,30,46)
    PROFILES_COR =bio_float.FloatSelector('CHLA', TI, R)
    wmo_list= bio_float.get_wmo_list(PROFILES_COR)
    wmo_list.sort()

    for wmo in wmo_list:
        WMO_PROFILES[wmo] = bio_float.filter_by_wmo(PROFILES_COR, wmo)


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
            WMO_PROFILES.setdefault(pCor._my_float.wmo, []).append(pCor)

# longest floats first, to balance the workers
wmo_order = sorted(WMO_PROFILES, key=lambda w: len(WMO_PROFILES[w]), reverse=True)
with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context('fork')) as executor:
    futures = {executor.submit(process_wmo, wmo): wmo for wmo in wmo_order}
    for future in as_completed(futures):
        wmo = futures[future]
        try:
            rows = future.result()
        except Exception as exc:
            print(f"{wmo} generated an exception: {exc!r}", flush=True)
            continue
        if rows:
            save_report(OUTDIAG, pd.concat(rows, ignore_index=True))
        print(f"{wmo} done: {len(rows)} profiles reported", flush=True)

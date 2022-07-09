import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates superfloat files of dissolved oxygen.
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
    parser.add_argument(   '--outdiag','-O',
                                type = str,
                                required = True,
                                default = "/gpfs/scratch/userexternal/gbolzon0/SUPERFLOAT/",
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
    return parser.parse_args()

args = argument()

if (args.datestart == 'NO_data') & (args.dateend == 'NO_data') & (args.update_file == 'NO_file'):
    raise ValueError("No file nor data inserted: you have to pass either datastart and dataeend or the update_file")

if ((args.datestart == 'NO_data') or (args.dateend == 'NO_data')) & (args.update_file == 'NO_file'):
    raise ValueError("No file nor data inserted: you have to pass both datastart and dataeend")


import pandas as pd
import CORIOLIS_checks
from instruments import bio_float
from Float import oxygen_saturation
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import superfloat_generator
from commons.utils import addsep
import os
import scipy.io.netcdf as NC
import numpy as np
import seawater as sw
from datetime import datetime, timedelta
from TREND_ANALYSIS import trend_conditions
from TREND_ANALYSIS import sign_analysis
import TREND_ANALYSIS
import basins.OGS as OGS
from instruments.var_conversions import FLOATVARS
from commons_local import cross_Med_basins, save_report

df_clim = pd.read_csv('EMODNET_climatology.csv',index_col=0)
df_cstd = pd.read_csv('EMODNET_stdev.csv',index_col=0)

class Metadata():
    def __init__(self, filename):
        self.filename = filename
        self.status_var = 'n'
        self.drift_code = -5
        self.offset = -999


def remove_bad_sensors(Profilelist,var):
    '''

    Subsetter, filtering out bad sensors for that var

     Arguments:
      * Profilelist * list of Profile objects
      * var         * string

      Returns:
        a list of Profile Objects
    '''
 
    OUT_N3n = ["6903197","6901767","6901773","6901771"]
    OUT_O2o = ["6901510"]
    OUT_O2o = ["6901766",'6903235','6902902',"6902700"]
    # 0 6901766 has negative values

    if ( var == 'SR_NO3' ):
        return [p for p in Profilelist if p.name() not in OUT_N3n]

    if ( var == 'DOXY' ):
        return [p for p in Profilelist if p.name() not in OUT_O2o]

    return Profilelist

def convert_oxygen(p,doxypres,doxyprofile):
    '''
    from micromol/Kg to  mmol/m3
    '''
    if doxypres.size == 0: return doxyprofile
    Pres, temp, Qc = p.read("TEMP",read_adjusted=False)
    Pres, sali, Qc = p.read("PSAL",read_adjusted=False)
    density = sw.dens(sali,temp,Pres)
    density_on_zdoxy = np.interp(doxypres,Pres,density)
    return doxyprofile * density_on_zdoxy/1000.

def dump_oxygen_file(outfile, p, Pres, Value, Qc, metadata, mode='w'):
    nP=len(Pres)
    if mode=='a':
        command = "cp %s %s.tmp" %(outfile,outfile)
        os.system(command)
    ncOUT = NC.netcdf_file(outfile + ".tmp" ,mode)

    if mode=='w': # if not existing file, we'll put header, TEMP, PSAL
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

    print("dumping oxygen on " + outfile, flush=True)
    doxy_already_existing="nDOXY" in ncOUT.dimensions.keys()
    if not doxy_already_existing : ncOUT.createDimension('nDOXY', nP)
    ncvar=ncOUT.createVariable("PRES_DOXY", 'f', ('nDOXY',))
    ncvar[:]=Pres
    ncvar=ncOUT.createVariable("DOXY", 'f', ('nDOXY',))
    ncvar[:]=Value
    #if not doxy_already_existing:
    setattr(ncvar, 'status_var' , metadata.status_var)
    setattr(ncvar, 'drift_code' , metadata.drift_code)
    setattr(ncvar, 'offset'     , metadata.offset)
    setattr(ncvar, 'variable'   , 'DOXY')
    setattr(ncvar, 'units'      , "mmol/m3")
    ncvar=ncOUT.createVariable("DOXY_QC", 'f', ('nDOXY',))
    ncvar[:]=Qc
    ncOUT.close()

    os.system("mv " + outfile + ".tmp " + outfile)


def get_outfile(p,outdir):
    wmo=p._my_float.wmo
    filename="%s%s/%s" %(outdir,wmo, os.path.basename(p._my_float.filename))
    return filename

def read_doxy(pCor):
    Pres, Value, Qc = pCor.read('DOXY',read_adjusted=True)
    nP=len(Pres)
    if nP<5 :
        print("few values for " + pCor._my_float.filename, flush=True)
        return None, None, None
    ValueCconv=convert_oxygen(pCor, Pres, Value)
    return Pres, ValueCconv, Qc

def Depth_interp(Profilelist, HistoricalDataset):
    """
    HistoricalDataset is a dictionary: keys are p.ID(), values are tuples (Pres,Value,Qc)

    data at 600m and 800m are created by interpolating data using layer 580-620m
        and 780-820 m
    """
    THRES             = 20
    LIST_DEPTH        = [600,800]
    COUNT=0
    VARNAME ,varmod   = 'DOXY' , 'O2o'
    columnlist        = ['ID','time','lat','lon','name','Type','VAR', 'Depth']
    df = pd.DataFrame(index=np.arange(0,2), columns=columnlist)
    for DEPTH in LIST_DEPTH:
        for profile in Profilelist:
            Pres, Profile, Qc = HistoricalDataset[profile.ID()] #read_doxy(profile)
            IDX = np.where((Pres >= DEPTH-THRES  ) & ( Pres <= DEPTH+THRES))
            lst = [profile.ID(), profile.time, profile.lat,profile.lon, profile.name(),
                  'Type', np.nan, np.nan]
            df.loc[COUNT] = pd.Series(lst , columnlist)
            if np.array(IDX).size ==0:
               df.Depth[COUNT] = DEPTH
            else:
               df.Depth[COUNT] = DEPTH
               df.VAR[COUNT]   = np.interp(DEPTH , Pres[IDX], (Profile[IDX]))
            COUNT+=1

    CONDITION = df[df.Depth==600].VAR.notnull().any()  # no data at 600m
    return(df, CONDITION)


def trend_analysis(p, Profilelist_hist, Dataset):
    starttime                  = p.time - timedelta(days=365*3)
    TI                         = TimeInterval.fromdatetimes(starttime, p.time)
    Profilelist                = [p for p in Profilelist_hist if TI.contains(p.time)]

    df, condition1_to_detrend  = Depth_interp(Profilelist, Dataset)
    ARGO       = Rectangle(np.float(p.lon) , np.float(p.lon) , np.float(p.lat) , np.float(p.lat))
    NAME_BASIN , BORDER_BASIN = cross_Med_basins(ARGO)
    return df, NAME_BASIN, condition1_to_detrend

def get_trend_report(p, df):
    LIST_DEPTH = [600,800]
    COLUMNS     = ['WMO','Depth','DURATION','min_date','max_date','Theil-Sen','RANSAC']
    df_report   = pd.DataFrame(index=np.arange(0,2), columns=COLUMNS)
    df_report['TREND_TIME_SERIES'] ,df_report['TREND_per_YEAR'] ,df_report['DRIFT_CODE'] = np.nan , np.nan , np.nan

    COUNT=0
    TI_3 = TimeInterval.fromdatetimes(p.time, p.time)

    for DEPTH in LIST_DEPTH:
        tmp = df[(df.Depth == DEPTH) & (df.name == wmo)]
        tmp.index = np.arange(0,len(tmp.index))
        days, min_d , max_d = CORIOLIS_checks.lenght_timeseries(tmp, 'time')
        Bool = CORIOLIS_checks.nans_check(tmp, 'VAR')
        tmp.dropna(inplace=True)
        lst = trend_conditions(wmo,days, Bool  , DEPTH,min_d, max_d , TI_3, tmp)
        df_report.iloc[COUNT,:] = pd.Series(lst, df_report.columns)
        serv = df_report.loc[df_report.WMO==wmo]
        A    = np.append(np.array( serv['Theil-Sen']), np.array( serv['RANSAC']))
        if serv.DURATION.any()==0:
            pass
        else:
            Bool = sign_analysis(A)
            df_report =  TREND_ANALYSIS.drift_coding(wmo, Bool, serv, df_report)
        COUNT+=1

    tmp = df[(df.Depth == 600) & (df.name == wmo)]
    return df_report, tmp

def clim_check(p, df_report, NAME_BASIN, tmp):
    VALCLIM    = float(df_clim.loc[df_clim.index==NAME_BASIN].iloc[:,0])
    TREND_null = df_report.TREND_TIME_SERIES.isnull().values.any()
    if TREND_null:
        OFFSET  = np.float(tmp.VAR.iloc[-1]) - VALCLIM
    else:
        Corrrected_val = np.float(tmp.VAR.iloc[-1]) - np.float(df_report.TREND_TIME_SERIES.iloc[0])
        OFFSET  = Corrrected_val - VALCLIM
    STDCLIM   = float(df_cstd.loc[df_cstd.index==NAME_BASIN].iloc[:,0])
    STDCLIM_2 = 2*STDCLIM
    wmo = p._my_float.wmo
    df_report.loc[ df_report.WMO == wmo , 'OFFSET'] = OFFSET
    df_report.loc[ df_report.WMO == wmo , 'basin'] = NAME_BASIN
    df_report.loc[ df_report.WMO == wmo , 'EMODNET'] = VALCLIM
    df_report.loc[ df_report.WMO == wmo , 'STdev*2'] =  STDCLIM_2
    threshold = STDCLIM_2
    return OFFSET, threshold, df_report

def apply_detrend(Pres, Prof_Coriolis, df_report):
    '''
    Arguments:
    * Pres          * ndarray
    * Prof_Coriolis * ndarray
    * df_report     *  dataframe

    Returns: Profile, ndarray - detrended oxygen

    Linear interp from 600m to 0
    Constant between 600m and bottom
    '''
    DEPTH=600
    Mask = [Pres<=DEPTH]
    Profile = Prof_Coriolis.copy()
    z_depth = np.array([Pres[0], DEPTH])
    x_var   = np.array([0, np.float(df_report.TREND_TIME_SERIES.iloc[0])])
    if np.isnan(np.sum(x_var)):
        pass
    else:
        fittval = np.polyfit(z_depth, x_var, 1)
        polyval = np.polyval(fittval, Pres[Mask]  )
        Profile[0:len(polyval)] = Profile[Mask]   -   polyval
        Mask_deep = [Pres>DEPTH]
        Profile[len(polyval):]  = Profile[Mask_deep] - polyval[-1]
    return Profile

def doxy_algorithm(p, Profilelist_hist, Dataset, outfile, metadata,writing_mode):
    '''
    Arguments:
    * p                * profile object
    * Profilelist_hist * List of profile object, provided by load_history()
    * Dataset          * dictionary, provided by load_history()
    '''
    Pres, Value, Qc  = Dataset[p.ID()]

#     Pres, Value, Qc = read_doxy(p)
#     if Pres is None: return
#
#
#
#     if p._my_float.status_var('DOXY')=='D':
#         metadata.status_var='D'
#         condition_to_write = True
#     else:
#         metadata.status_var='A'
#         condition_to_write =  oxygen_saturation.oxy_check(Pres,Value,p)
#
#
#     if not condition_to_write:
#         print("Saturation Test not passed")
#         return
    metadata.status_var = p._my_float.status_var('DOXY')
    df, NAME_BASIN, condition1_to_detrend = trend_analysis(p, Profilelist_hist, Dataset)
    Oxy_Profile = Value

    if not condition1_to_detrend:
        DRIFT_CODE = -1
        print("no detrend possible", flush=True)
        metadata.drift_code = DRIFT_CODE
    else:
        df_report, tmp = get_trend_report(p, df)
        if tmp[tmp.Depth==600].VAR.isnull().values.any() :
            DRIFT_CODE = -1
            metadata.drift_code = DRIFT_CODE
        else:

            OFFSET, threshold, df_report = clim_check(p,df_report, NAME_BASIN, tmp)

            if abs(OFFSET) >= threshold:
                wmo = p._my_float.wmo
                df_report.loc[ (df_report.WMO == wmo) & (df_report.Depth== 600), 'Black_list'] = 'True'
                timenum = int(p.time.strftime("%Y%m%d"))
                save_report( OUT_META+ "Blacklist_wmo.csv", 1,['WMO', 'DATE_DAY' , 'OFFSET' , 'STDCLIM_2'],[int(wmo), timenum, OFFSET , threshold])
                return
            else:
                Oxy_Profile = apply_detrend(Pres, Value, df_report)
                metadata.drift_code = df_report['DRIFT_CODE'].iloc[0]
                metadata.offset = OFFSET
                # high freq.csv

    os.system('mkdir -p ' + os.path.dirname(outfile))
    dump_oxygen_file(outfile, p, Pres, Oxy_Profile, Qc, metadata,mode=writing_mode)



def load_history(wmo):
    '''
     Replicates superfloat dataset without detrend - the previous doxy_algorithm
    '''
    print("Loading dataset for float", wmo, "...", flush=True)
    TI     = TimeInterval("1950","2050",'%Y')
    R = Rectangle(-6,36,30,46)
    PROFILES_COR_all =bio_float.FloatSelector('DOXY', TI, R)
    PROFILES_COR = remove_bad_sensors(PROFILES_COR_all, "DOXY")
    Profilelist_wmo=bio_float.filter_by_wmo(PROFILES_COR, wmo)
    Profilelist=[]
    Dataset={}
    for p in Profilelist_wmo:
        if p._my_float.status_var('DOXY')=='R': continue
        Pres, Value, Qc = read_doxy(p)
        if Pres is None: continue
        if p._my_float.status_var('DOXY')=='D':
            condition_to_write = True
        else:
            condition_to_write =  oxygen_saturation.oxy_check(Pres,Value,p)
        if not condition_to_write:
            print(p._my_float.filename, "Saturation Test not passed", flush=True)
            continue

        Profilelist.append(p)
        Dataset[p.ID()] = (Pres,Value,Qc)
    print("...done", flush=True)
    return Profilelist, Dataset



OUTDIR = addsep(args.outdir)
OUT_META = addsep(args.outdiag)
input_file=args.update_file
if input_file == 'NO_file':
    TI     = TimeInterval(args.datestart,args.dateend,'%Y%m%d')
    R = Rectangle(-6,36,30,46)

    PROFILES_COR_all =bio_float.FloatSelector('DOXY', TI, R)
    PROFILES_COR = remove_bad_sensors(PROFILES_COR_all, "DOXY")

    wmo_list= bio_float.get_wmo_list(PROFILES_COR)
    wmo_list.sort()


    for wmo in wmo_list:
        print (wmo, flush=True)

        Hist_filtered_Profilelist, Dataset = load_history(wmo)
        Selected_Profilelist=bio_float.filter_by_wmo(PROFILES_COR, wmo)
        Profilelist= [p for p in Selected_Profilelist if p in Hist_filtered_Profilelist]
        for ip, p in enumerate(Profilelist):
            outfile = get_outfile(p,OUTDIR)

            writing_mode=superfloat_generator.writing_mode(outfile)

            condition_to_write = ~superfloat_generator.exist_valid_variable('DOXY',outfile)
            if args.force: condition_to_write=True
            if not condition_to_write: continue

            metadata = Metadata(p._my_float.filename)
            doxy_algorithm(p, Hist_filtered_Profilelist, Dataset , outfile, metadata,writing_mode)
else:

    INDEX_FILE=superfloat_generator.read_float_update(input_file)
    nFiles=INDEX_FILE.size

    for iFile in range(nFiles):
        timestr          = INDEX_FILE['date'][iFile].decode()
        lon              = INDEX_FILE['longitude' ][iFile]
        lat              = INDEX_FILE['latitude' ][iFile]
        filename         = INDEX_FILE['file_name'][iFile].decode()
        available_params = INDEX_FILE['parameters'][iFile].decode()
        parameterdatamode= INDEX_FILE['parameter_data_mode'][iFile].decode()
        float_time = datetime.strptime(timestr,'%Y%m%d%H%M%S')
        filename=filename.replace('coriolis/','').replace('profiles/','')


        if 'DOXY' not in available_params: continue
        p=bio_float.profile_gen(lon, lat, float_time, filename, available_params,parameterdatamode)
        wmo=p._my_float.wmo
        OUT_O2o = ["6901766",'6903235','6902902',"6902700"]

        if wmo in OUT_O2o: continue
        outfile = get_outfile(p,OUTDIR)
        if p._my_float.status_var('DOXY')=='R': continue

        writing_mode=superfloat_generator.writing_mode(outfile)
        metadata = Metadata(p._my_float.filename)
        doxy_algorithm(p, outfile, metadata, writing_mode)

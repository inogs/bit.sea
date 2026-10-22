import argparse
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from bitsea.utilities.argparse_types import existing_dir_path,existing_file_path
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
                                type = existing_dir_path,
                                required = True,
                                default = "/gpfs/scratch/userexternal/gbolzon0/SUPERFLOAT/",
                                help = 'path of the Superfloat dataset ')
    parser.add_argument(   '--outdiag','-O',
                                type = existing_dir_path,
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
    parser.add_argument(   '--make_plots',
                                action = argparse.BooleanOptionalAction,
                                default = True,
                                help = '''Generate, for each float where a drift correction was applied,
                                a diagnostic figure comparing the raw O2 timeseries (red) with the
                                drift-corrected QC_O2 timeseries (blue) at 600m and 800m.
                                Enabled by default; disable with --no-make_plots.''')
    parser.add_argument(   '--tmp_report_dir','-T',
                                type = str,
                                required = False,
                                default = None,
                                help = 'path for temporary worker report files')
    parser.add_argument(   '--workers','-w',
                                type = int,
                                required = False,
                                default = os.cpu_count() or 1,
                                help = 'number of worker processes for parallel execution')
    return parser.parse_args()


args = argument()

if (args.datestart == 'NO_data') & (args.dateend == 'NO_data') & (args.update_file == 'NO_file'):
    raise ValueError("No file nor data inserted: you have to pass either datestart and dataeend or the update_file")

if ((args.datestart == 'NO_data') or (args.dateend == 'NO_data')) & (args.update_file == 'NO_file'):
    raise ValueError("No file nor data inserted: you have to pass both datestart and dateend or update_file")


import pandas as pd
import CORIOLIS_checks
from bitsea.instruments import bio_float
from bitsea.Float import oxygen_saturation
from bitsea.commons.time_interval import TimeInterval
from bitsea.basins.region import Rectangle
import superfloat_generator
from pathlib import Path
import os
import netCDF4 as NC
import numpy as np
import gsw
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import TREND_ANALYSIS
import bitsea.basins.OGS as OGS
from bitsea.instruments.var_conversions import FLOATVARS
from commons_local import cross_Med_basins

df_clim = pd.read_csv('EMODNET_climatology.csv',index_col=0)
df_cstd = pd.read_csv('EMODNET_stdev.csv',index_col=0)

# Annual drift-rate thresholds (mmol/m3/yr), computed once per float on its
# entire historical timeseries (see compute_float_trend()):
#  * |TREND_per_YEAR| <  TREND_THRESHOLD_MIN -> not significant, no correction
#  * |TREND_per_YEAR| >  TREND_THRESHOLD_MAX -> anomalous drift, float rejected
#  * otherwise                               -> drift correction applied
TREND_THRESHOLD_MIN = 2.0
TREND_THRESHOLD_MAX = 10.0

class Metadata():
    def __init__(self, filename):
        self.filename = str(filename)
        self.status_var = 'n'
        self.drift_code = np.nan
        self.offset = np.nan
        self.drift  = np.nan


def read_temp_psal(p):
    PresT, Temp, QcT = p.read('TEMP', read_adjusted=True)
    Pres, Sali, QcS = p.read('PSAL', read_adjusted=True)
    if (Pres is None or PresT is None or Temp is None or Sali is None or len(Pres) < 5 or len(PresT) < 5):
        PresT, Temp, QcT = p.read('TEMP', read_adjusted=False)
        Pres, Sali, QcS = p.read('PSAL', read_adjusted=False)
    return PresT, Temp, QcT, Pres, Sali, QcS


def convert_oxygen(p,doxypres,doxyprofile):
    '''
    from micromol/Kg to  mmol/m3
    '''
    if doxypres.size == 0: return doxyprofile
    PresT, temp, Qc, Pres, sali, QcS = read_temp_psal(p)
    if len(temp) != len(sali):
           temp = np.interp(Pres, PresT, temp)
    SA = gsw.SA_from_SP(sali, Pres, p.lon, p.lat)
    density = gsw.rho(SA, gsw.CT_from_t(SA, temp, Pres), Pres)
    density_on_zdoxy = np.interp(doxypres,Pres,density)
    return doxyprofile * density_on_zdoxy/1000.

def save_df_report(out_dir, df_report, filename, tmp_report_dir=None, wmo=None):
    """
    Save or update a CSV report file by appending the first row of a DataFrame.

    Supports temporary per-worker report output when tmp_report_dir is provided.
    """

    if tmp_report_dir is not None:
        out_dir = Path(tmp_report_dir)
        if wmo is not None:
            filename = f"{filename}.{wmo}.csv"
    else:
        out_dir = Path(out_dir)

    out_dir.mkdir(parents=True, exist_ok=True)
    file_path = out_dir / filename

    df_new = df_report.iloc[[0]].copy()
    if file_path.exists():
        df_existing = pd.read_csv(file_path)
        df_combined = pd.concat([df_existing, df_new], ignore_index=True)
    else:
        df_combined = df_new

    df_combined = df_combined.drop_duplicates()
    df_combined.to_csv(file_path, index=False)


def combine_temp_reports(tmp_report_dir, out_meta, filename):
    tmp_report_dir = Path(tmp_report_dir)
    out_meta = Path(out_meta)
    out_meta.mkdir(parents=True, exist_ok=True)
    output_file = out_meta / filename

    temp_files = sorted(tmp_report_dir.glob(f"{filename}.*.csv"))
    if not temp_files:
        return

    dfs = [pd.read_csv(f) for f in temp_files]
    if output_file.exists():
        dfs.insert(0, pd.read_csv(output_file))
    df_combined = pd.concat(dfs, ignore_index=True).drop_duplicates()
    df_combined.to_csv(output_file, index=False)

def dump_oxygen_file(outfile, p, Pres, Value, Qc, metadata, mode='w'):
    nP=len(Pres)
    tmpfile=outfile.with_suffix('.nc.tmp')
    if mode=='a':
        command = "cp {} {}".format(outfile,tmpfile)
        os.system(command)
    with NC.Dataset(tmpfile,mode) as ncOUT:

        if mode=='w': # if not existing file, we'll put header, TEMP, PSAL
            setattr(ncOUT, 'origin'     , 'coriolis')
            setattr(ncOUT, 'file_origin', metadata.filename)
            PresT, Temp, QcT, PresS, Sali, QcS = read_temp_psal(p)
            ncOUT.createDimension("DATETIME",14)
            ncOUT.createDimension("NPROF", 1)
            ncOUT.createDimension('nTEMP', len(PresT))
            ncOUT.createDimension('nPSAL', len(PresS))

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

            ncvar=ncOUT.createVariable('PSAL','f',('nPSAL',))
            ncvar[:]=Sali
            setattr(ncvar, 'variable'   , 'SALI')
            setattr(ncvar, 'units'      , "PSS78")
            ncvar=ncOUT.createVariable('PRES_PSAL','f',('nPSAL',))
            ncvar[:]=PresS
            ncvar=ncOUT.createVariable('PSAL_QC','f',('nPSAL',))
            ncvar[:]=QcS

        print("dumping oxygen on " + str(outfile), flush=True)
        doxy_already_existing=superfloat_generator.exist_valid_variable("DOXY", tmpfile)
        if not doxy_already_existing :
            ncOUT.createDimension('nDOXY', nP)
            ncvar=ncOUT.createVariable("PRES_DOXY", 'f', ('nDOXY',))
            ncvar[:]=Pres
            ncvar=ncOUT.createVariable("DOXY", 'f', ('nDOXY',))
            ncvar[:]=Value
            ncvar=ncOUT.createVariable("DOXY_QC", 'f', ('nDOXY',))
            ncvar[:]=Qc
        else:
            ncvar=ncOUT.variables['PRES_DOXY']
            ncvar[:]=Pres
            ncvar=ncOUT.variables['DOXY']
            ncvar[:]=Value
            ncvar=ncOUT.variables['DOXY_QC']
            ncvar[:]=Qc
        setattr(ncvar, 'status_var' , metadata.status_var)
        setattr(ncvar, 'drift_code' , metadata.drift_code)
        setattr(ncvar, 'drift'      , metadata.drift)
        setattr(ncvar, 'offset'     , metadata.offset)
        setattr(ncvar, 'variable'   , 'DOXY')
        setattr(ncvar, 'units'      , "mmol/m3")

    os.system("mv {} {}".format(tmpfile,outfile))


def get_outfile(p,outdir):
    wmo=p._my_float.wmo
    filename = outdir / wmo / p._my_float.filename.name
    return filename

def read_doxy(pCor):
    Pres, Value, Qc = pCor.read('DOXY',read_adjusted=True)
    nP=len(Pres)
    if nP<5 :
        print("few values for " + str(pCor._my_float.filename), flush=True)
        return None, None, None
    ValueCconv=convert_oxygen(pCor, Pres, Value)
    return Pres, ValueCconv, Qc

def Depth_interp(Profilelist, HistoricalDataset):
    """
    Interpolate oxygen (DOXY) values at standard depths (600m and 800m) from historical profiles.

    Parameters
    ----------
    Profilelist : list of Profile profiles to process.
    HistoricalDataset : dict where keys are profile IDs and values are tuples (Pres, Value, Qc) containing pressure, measured oxygen values, and quality control flags.

    What it does
    ------------
    1. For each profile in Profilelist, checks if data exists within ±20m around 600m and 800m.
    2. Interpolates oxygen (DOXY) values at exactly 600m and 800m if sufficient data exists.
    3. If the current profile lacks valid data at a target depth, it may reuse the most recent
       valid value from up to two previous profiles within a short temporal window.
    4. Creates a DataFrame summarizing profile ID, time, location, depth, and interpolated values.
    5. Sets a condition flag True if at least one valid oxygen value exists at 600m.

    Returns
    -------
    df : DataFrame with interpolated oxygen values at 600m and 800m for each profile.
    CONDITION : bool True if at least one valid oxygen value exists at 600m, False otherwise.
    """

    THRES             = 20  # data are interpolated using layer +/- 20m around target depth
    FALLBACK_MAX_DAYS = 10
    MAX_PREVIOUS      = 2
    LIST_DEPTH        = [600, 800]
    columnlist        = ['ID', 'time', 'lat', 'lon', 'name', 'Type', 'VAR', 'Depth']
    df = pd.DataFrame(index=np.arange(0, 2 * len(Profilelist)), columns=columnlist)
    Profilelist_sorted = sorted(Profilelist, key=lambda p: p.time)
    COUNT = 0

    for DEPTH in LIST_DEPTH:
        previous_valid = []  # list of tuples (profile_time, interpolated_value)
        for profile in Profilelist_sorted:
            Pres, Profile, Qc = HistoricalDataset[profile.ID()]
            IDX = np.where((Pres >= DEPTH-THRES) & (Pres <= DEPTH+THRES))
            lst = [profile.ID(), profile.time, profile.lat, profile.lon, profile.name(),
                   'Type', np.nan, DEPTH]
            df.loc[COUNT] = pd.Series(lst, columnlist)

            value = np.nan
            if IDX[0].size != 0:
                value = np.interp(DEPTH, Pres[IDX], Profile[IDX])
            else:
                for prev_time, prev_value in reversed(previous_valid[-MAX_PREVIOUS:]):
                    day_diff = (profile.time - prev_time).days
                    if day_diff <= FALLBACK_MAX_DAYS:
                        value = prev_value
                        break
                    if day_diff > FALLBACK_MAX_DAYS:
                        break

            if not np.isnan(value):
                df.at[COUNT, 'VAR'] = value
                previous_valid.append((profile.time, value))

            COUNT += 1

    CONDITION = df[df.Depth == 600].VAR.notnull().any()
    return df, CONDITION


def get_basin(p):
    '''
    Returns the name of the Mediterranean sub-basin containing profile p.
    '''
    ARGO       = Rectangle(float(p.lon) , float(p.lon) , float(p.lat) , float(p.lat))
    NAME_BASIN , BORDER_BASIN = cross_Med_basins(ARGO)
    return NAME_BASIN

def compute_float_trend(wmo, Profilelist_hist, Dataset):
    """
    Computes a single dissolved-oxygen drift rate for the WHOLE record of a
    BGC-Argo float (not a per-profile rolling window).

    Parameters
    ----------
    wmo              : str, WMO code of the float
    Profilelist_hist : list of Profile, entire (QC-filtered) history of the float
    Dataset          : dict, {profile.ID(): (Pres, Value, Qc)}, from load_history()

    What it does
    ------------
    1. Interpolates DOXY at 600m and 800m (+/-20m) for every profile in the
       float's history (Depth_interp).
    2. Requires the series to span >= 365 days and to have valid 600m data,
       otherwise no trend is computed (DRIFT_CODE=-1).
    3. Fits Theil-Sen and RANSAC regressors at each depth on the REAL elapsed
       time (days since the float's first profile), not on the sample index.
    4. Checks sign consistency of the 4 estimates (sign_analysis).
    5. Classifies the float:
         abs(TREND_per_YEAR) below TREND_THRESHOLD_MIN -> DRIFT_CODE = 0 (not significant)
         abs(TREND_per_YEAR) above TREND_THRESHOLD_MAX -> DRIFT_CODE = -2 (anomalous, rejected)
         otherwise                                     -> DRIFT_CODE = 1 (correction applied)

    Returns
    -------
    float_trend : dict, one row of information at float level (see keys below)
    df_hist     : pandas.DataFrame, 600m/800m interpolated series for the whole
                  float history (reused for the climatology check and the plots)
    """

    df_hist, condition1 = Depth_interp(Profilelist_hist, Dataset)

    t0   = min(p.time for p in Profilelist_hist)
    tmax = max(p.time for p in Profilelist_hist)
    duration_days = (tmax - t0).days

    float_trend = {
        'WMO': wmo, 't0': t0, 'series_start': t0, 'series_end': tmax,
        'duration_days': duration_days,
        'TheilSen_600': np.nan, 'RANSAC_600': np.nan,
        'TheilSen_800': np.nan, 'RANSAC_800': np.nan,
        'sign_consistent': False, 'TREND_per_YEAR': np.nan, 'DRIFT_CODE': -1,
    }

    if (not condition1) or (duration_days < 365):
        return float_trend, df_hist

    trend_by_depth = {}
    for DEPTH in [600, 800]:
        tmp = df_hist[(df_hist.Depth == DEPTH) & (df_hist.name == wmo)].dropna(subset=['VAR']).copy()
        if len(tmp) < 5:
            trend_by_depth[DEPTH] = (np.nan, np.nan)
            continue
        tmp['time'] = pd.to_datetime(tmp['time'])
        x = np.array([(t - t0).days for t in tmp['time']], dtype=float)
        y = np.array(tmp['VAR'], dtype=float)
        estimators = TREND_ANALYSIS.Trend_RANSAC_ThSEN(x, y)
        trend_by_depth[DEPTH] = (estimators[0][1], estimators[1][1])

    float_trend['TheilSen_600'], float_trend['RANSAC_600'] = trend_by_depth[600]
    float_trend['TheilSen_800'], float_trend['RANSAC_800'] = trend_by_depth[800]

    A = np.array([float_trend['TheilSen_600'], float_trend['TheilSen_800'],
                  float_trend['RANSAC_600'], float_trend['RANSAC_800']], dtype=float)

    if np.isnan(A).any() or not TREND_ANALYSIS.sign_analysis(A):
        return float_trend, df_hist

    float_trend['sign_consistent'] = True
    AVG_TREND = np.mean(A)
    YR_TREND  = (AVG_TREND * 365) / duration_days
    float_trend['TREND_per_YEAR'] = YR_TREND

    if abs(YR_TREND) < TREND_THRESHOLD_MIN:
        float_trend['DRIFT_CODE'] = 0
    elif abs(YR_TREND) > TREND_THRESHOLD_MAX:
        float_trend['DRIFT_CODE'] = -2
    else:
        float_trend['DRIFT_CODE'] = 1

    return float_trend, df_hist

def clim_check(p, df_hist, wmo, NAME_BASIN, correction_600m, yr_trend):
    """
    Compare, for a single profile, the most recent 600m oxygen value (up to
    and including p.time) against the EMODNET climatology of its basin.

    If the float has a significant drift trend, the observed value is
    corrected (subtracting correction_600m, the progressive drift correction
    already computed for this profile) before the comparison.

    Parameters
    ----------
    p               : Profile, current profile.
    df_hist         : pandas.DataFrame, 600m/800m interpolated series for the whole float history.
    wmo             : str, WMO code of the float.
    NAME_BASIN      : str, Mediterranean sub-basin of this profile.
    correction_600m : float, drift correction (elapsed_years * yr_trend) for this profile.
    yr_trend        : float, float-level annual drift rate (TREND_per_YEAR).

    Returns
    -------
    OFFSET    : float, observed (or corrected) value minus EMODNET climatology.
    threshold : float, climatological threshold (2 x stdev, floored at 15).
    VALCLIM   : float, EMODNET climatology used for this basin.
    STDCLIM_2 : float, threshold actually used (== threshold, kept for reporting).
    """
    VALCLIM = float(df_clim.loc[df_clim.index == NAME_BASIN].iloc[:, 0])

    tmp600 = df_hist[(df_hist.Depth == 600) & (df_hist.name == wmo)].copy()
    tmp600['time'] = pd.to_datetime(tmp600['time'])
    tmp600 = tmp600[tmp600['time'] <= p.time]
    vals = tmp600['VAR'].dropna().tail(3)
    if vals.empty:
        raise ValueError("clim_check(): no valid 600m values available for OFFSET computation")

    observed = float(vals.iloc[-1])
    if np.isnan(yr_trend) or abs(yr_trend) < TREND_THRESHOLD_MIN:
        OFFSET = observed - VALCLIM
    else:
        OFFSET = (observed - correction_600m) - VALCLIM

    STDCLIM   = float(df_cstd.loc[df_cstd.index == NAME_BASIN].iloc[:, 0])
    STDCLIM_2 = 2 * STDCLIM
    if STDCLIM_2 < 15:
        STDCLIM_2 = 15
    threshold = STDCLIM_2
    return OFFSET, threshold, VALCLIM, STDCLIM_2

def apply_detrend(Pres, Prof_Coriolis, drift_600m):
    '''
    Arguments:
    * Pres          * ndarray
    * Prof_Coriolis * ndarray
    * drift_600m    * float, drift accumulated up to this profile's date
                      (elapsed_years * TREND_per_YEAR), applied at 600m

    Returns: Profile, ndarray - detrended oxygen

    Linear interp from 600m to 0
    Constant between 600m and bottom
    '''
    DEPTH = 600
    Mask = Pres <= DEPTH
    Profile = Prof_Coriolis.copy()
    if np.isnan(drift_600m):
        return Profile
    if not np.any(Mask):
        return Profile
    z_depth = np.array([Pres[0], DEPTH])
    x_var = np.array([0, drift_600m])
    fittval = np.polyfit(z_depth, x_var, 1)
    polyval = np.polyval(fittval, Pres[Mask])
    if len(polyval) == 0:
        return Profile
    Profile[0:len(polyval)] = Profile[Mask] - polyval
    Mask_deep = Pres > DEPTH
    Profile[len(polyval):] = Profile[Mask_deep] - polyval[-1]
    return Profile


def plot_drift_correction(wmo, df_hist, float_trend, outdir):
    '''
    Saves one figure per float (2 rows x 1 column: 600m on top, 800m below)
    comparing the raw O2 timeseries (red) with the drift-corrected QC_O2
    timeseries (blue), for the whole float history.
    Only meaningful when float_trend['DRIFT_CODE'] == 1.
    '''
    t0       = float_trend['t0']
    yr_trend = float_trend['TREND_per_YEAR']
    plot_dir = Path(outdir) / "PLOTS_DRIFT"
    plot_dir.mkdir(parents=True, exist_ok=True)

    fig, axs = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    for ax, DEPTH in zip(axs, [600, 800]):
        tmp = df_hist[(df_hist.Depth == DEPTH) & (df_hist.name == wmo)].dropna(subset=['VAR']).copy()
        if tmp.empty:
            ax.set_title(f"{DEPTH} m - no data")
            continue
        tmp['time'] = pd.to_datetime(tmp['time'])
        tmp = tmp.sort_values('time')
        elapsed_years = (tmp['time'] - t0).dt.days / 365.25
        capped_years = np.minimum(elapsed_years, 1.0)
        corrected = tmp['VAR'] - yr_trend * capped_years
        ax.plot(tmp['time'], tmp['VAR'], 'r', label='O2')
        ax.plot(tmp['time'], corrected, 'b', label='QC_O2')
        ax.set_ylabel('O2 (mmol/m3)')
        ax.set_title(f'{DEPTH} m')
        ax.legend(loc='upper left')
        ax.grid(True)
    axs[-1].tick_params(labelrotation=45)
    fig.suptitle(f'WMO {wmo} - drift correction (trend = {yr_trend:.2f} mmol/m3/yr)')
    fig.tight_layout()
    fig.savefig(plot_dir / f"{wmo}.png")
    plt.close(fig)


def build_trend_report_row(float_trend):
    return pd.DataFrame([{
        'WMO'            : float_trend['WMO'],
        'series_start'   : float_trend['series_start'],
        'series_end'     : float_trend['series_end'],
        'duration_days'  : float_trend['duration_days'],
        'TheilSen_600'   : float_trend['TheilSen_600'],
        'RANSAC_600'     : float_trend['RANSAC_600'],
        'TheilSen_800'   : float_trend['TheilSen_800'],
        'RANSAC_800'     : float_trend['RANSAC_800'],
        'sign_consistent': float_trend['sign_consistent'],
        'TREND_per_YEAR' : float_trend['TREND_per_YEAR'],
        'DRIFT_CODE'     : float_trend['DRIFT_CODE'],
    }])


def build_accepted_row(wmo, p, NAME_BASIN, float_trend, correction_600m, OFFSET, VALCLIM, STDCLIM_2):
    return pd.DataFrame([{
        'WMO'                    : wmo,
        'cycle'                  : p._my_float.cycle,
        'time'                   : p.time,
        'basin'                  : NAME_BASIN,
        'DRIFT_CODE'             : float_trend['DRIFT_CODE'],
        'TREND_per_YEAR'         : float_trend['TREND_per_YEAR'],
        'elapsed_years'          : (p.time - float_trend['t0']).days / 365.25,
        'correction_600m_applied': correction_600m,
        'OFFSET'                 : OFFSET,
        'EMODNET'                : VALCLIM,
        'STdev*2'                : STDCLIM_2,
    }])


def build_rejected_row(wmo, p, NAME_BASIN, OFFSET, VALCLIM, STDCLIM_2, TREND_per_YEAR, reject_reason):
    return pd.DataFrame([{
        'WMO'           : wmo,
        'cycle'         : p._my_float.cycle,
        'time'          : p.time,
        'basin'         : NAME_BASIN,
        'reject_reason' : reject_reason,
        'OFFSET'        : OFFSET,
        'EMODNET'       : VALCLIM,
        'STdev*2'       : STDCLIM_2,
        'TREND_per_YEAR': TREND_per_YEAR,
    }])

def doxy_algorithm(p, wmo, float_trend, df_hist, Dataset, outfile, metadata, writing_mode, out_meta, tmp_report_dir=None):
    '''
    Arguments:
    * p           * profile object (current profile to write)
    * wmo         * str, WMO code of the float
    * float_trend * dict, float-level trend info from compute_float_trend()
    * df_hist     * DataFrame, 600m/800m interpolated series for the whole float history
    * Dataset     * dictionary, provided by load_history()
    * out_meta    * path where metadata reports are written
    * tmp_report_dir * optional path for temporary worker report files
    '''
    Pres, Value, Qc  = Dataset[p.ID()]
    PresT, _, _, _, _, _ = read_temp_psal(p)
    if PresT is None or len(PresT) < 5:
        print("few values in Coriolis TEMP in " + str(p._my_float.filename), flush=True)
        return

    metadata.status_var = p._my_float.status_var('DOXY')
    NAME_BASIN = get_basin(p)
    Oxy_Profile = Value

    DRIFT_CODE = float_trend['DRIFT_CODE']
    YR_TREND   = float_trend['TREND_per_YEAR']
    correction_600m = np.nan
    if DRIFT_CODE == 1:
        elapsed_years = (p.time - float_trend['t0']).days / 365.25
        capped_years = min(elapsed_years, 1.0)
        correction_600m = YR_TREND * capped_years

    no_clim_value = False
    try:
        OFFSET, threshold, VALCLIM, STDCLIM_2 = clim_check(p, df_hist, wmo, NAME_BASIN, correction_600m, YR_TREND)
    except ValueError:
        print("no valid 600m climatology reference for " + str(p._my_float.filename), flush=True)
        OFFSET = np.nan
        VALCLIM = np.nan
        STDCLIM_2 = np.nan
        no_clim_value = True

    if not no_clim_value and abs(OFFSET) >= threshold: # discard profile
        row = build_rejected_row(wmo, p, NAME_BASIN, OFFSET, VALCLIM, STDCLIM_2, YR_TREND, 'climatology_offset')
        save_df_report(out_meta, row, "Floats_rejected.csv", tmp_report_dir=tmp_report_dir, wmo=wmo)
        return

    if no_clim_value:
        metadata.drift_code = -1
        metadata.offset = np.nan
        metadata.drift = np.nan
    elif DRIFT_CODE == 1 and np.any(Pres <= 600):
        Oxy_Profile = apply_detrend(Pres, Value, correction_600m)
        metadata.drift_code = 1
        metadata.offset = OFFSET
        metadata.drift  = YR_TREND
    else:
        metadata.drift_code = DRIFT_CODE
        metadata.offset = OFFSET

    row = build_accepted_row(wmo, p, NAME_BASIN, float_trend, correction_600m, OFFSET, VALCLIM, STDCLIM_2)
    save_df_report(out_meta, row, "Floats_accepted.csv", tmp_report_dir=tmp_report_dir, wmo=wmo)
    Path.mkdir(outfile.parent, exist_ok=True)
    dump_oxygen_file(outfile, p, Pres, Oxy_Profile, Qc, metadata,mode=writing_mode)



def load_history(wmo, write_report=False, out_meta=None, tmp_report_dir=None):
    '''
     Replicates superfloat dataset without detrend - the previous doxy_algorithm
     Adds optional logging to a CSV in OUT_META or temporary report directory.

     Input: wmo          : str , WMO code of the float to load
            write_report : bool, If True, rejected profiles (status 'R', missing data, saturation test fail) are logged.
            out_meta     : path where final metadata reports are saved.
            tmp_report_dir : optional path for temporary worker report files.

     Output: Profilelist : List of profiles that passed Data Mode and Saturation checks
             Dataset     : Dictionary of data keyed by profile ID, values = (Pres, Value, Qc)
    '''

    df_reject = None
    if write_report:
        if tmp_report_dir is not None:
            tmp_report_dir = Path(tmp_report_dir)
            tmp_report_dir.mkdir(parents=True, exist_ok=True)
            csv_file = tmp_report_dir / f"DataMode_and_Saturation_rejection_doxy.csv.{wmo}.csv"
            df_reject = pd.DataFrame(columns=["wmo","cycle","date","reasons","basin"])
        else:
            out_meta = Path(out_meta) if out_meta is not None else Path('.')
            out_meta.mkdir(parents=True, exist_ok=True)
            csv_file = out_meta / "DataMode_and_Saturation_rejection_doxy.csv"
            if csv_file.exists():
                df_reject = pd.read_csv(csv_file)
            else:
                df_reject = pd.DataFrame(columns=["wmo","cycle","date","reasons","basin"])

        def log_reject(wmo_name, cycle,date, reason, NAME_BASIN):
            nonlocal df_reject
            df_reject = pd.concat(
                [df_reject, pd.DataFrame({"wmo": [wmo_name], "cycle": [cycle] , "date": [date], "reasons": [reason], "basin": [NAME_BASIN] })],
                ignore_index=True)
            df_reject.drop_duplicates(inplace=True)

    print("Loading dataset for float", wmo, "...", flush=True)

    TI     = TimeInterval("1950","2050",'%Y')
    R = Rectangle(-6,36,30,46)
    PROFILES_COR =bio_float.FloatSelector('DOXY', TI, R)
    Profilelist_wmo=bio_float.filter_by_wmo(PROFILES_COR, wmo)
    Profilelist=[]
    Dataset={}
    for p in Profilelist_wmo:
        
        ARGO       = Rectangle(float(p.lon) , float(p.lon) , float(p.lat) , float(p.lat))
        NAME_BASIN , BORDER_BASIN = cross_Med_basins(ARGO)
        status = p._my_float.status_var('DOXY')
        
        if status == 'R':
            if write_report:
                log_reject(p.name(), p._my_float.cycle, p.time, "RT" ,NAME_BASIN )
            continue
        
        Pres, Value, Qc = read_doxy(p) # convert doxy and filter len pres<5
        if Pres is None:
            if write_report:
                log_reject(p.name(), p._my_float.cycle, p.time, "PresNone", NAME_BASIN)
            continue
        
        if status == 'D':
            condition_to_write = True
        else: #ADJUSTED MODE
            condition_to_write = oxygen_saturation.oxy_check(Pres, Value, p)
        
        if not condition_to_write:
            if write_report:
                log_reject(p.name(), p._my_float.cycle, p.time, "SaturationTest", NAME_BASIN)
            continue
        Profilelist.append(p)
        Dataset[p.ID()] = (Pres,Value,Qc)

    if write_report and df_reject is not None:
        df_reject.to_csv(csv_file, index=False)
    print("...done", flush=True)
    return Profilelist, Dataset


OUTDIR = Path(args.outdir)
OUT_META = Path(args.outdiag)


def process_wmo(wmo, outdir, out_meta, args_dict):
    outdir = Path(outdir)
    out_meta = Path(out_meta)
    tmp_report_dir = Path(args_dict['tmp_report_dir']) if args_dict.get('tmp_report_dir') else out_meta / 'tmp_reports'
    tmp_report_dir.mkdir(parents=True, exist_ok=True)

    TI = TimeInterval(args_dict['datestart'], args_dict['dateend'], '%Y%m%d')
    R = Rectangle(-6,36,30,46)
    PROFILES_COR_local = bio_float.FloatSelector('DOXY', TI, R)

    Hist_filtered_Profilelist, Dataset = load_history(wmo, True, out_meta=out_meta, tmp_report_dir=tmp_report_dir)
    Selected_Profilelist = bio_float.filter_by_wmo(PROFILES_COR_local, wmo)
    Profilelist = [p for p in Selected_Profilelist if p in Hist_filtered_Profilelist]
    if not Profilelist:
        return {'wmo': wmo, 'skipped': True, 'message': 'no valid profiles'}

    profile_doxy_exists = {}
    profile_file_valid = {}
    for p in Profilelist:
        outfile = get_outfile(p, outdir)
        profile_file_valid[p.ID()] = superfloat_generator.exist_valid(outfile)
        profile_doxy_exists[p.ID()] = profile_file_valid[p.ID()] and superfloat_generator.exist_variable('DOXY', outfile)

    existing_count = sum(profile_doxy_exists.values())
    total_count = len(Profilelist)
    print(f"{wmo}: {existing_count}/{total_count} profili già presenti in {outdir}", flush=True)
    if existing_count == total_count:
        print(f"{wmo} già completo in {outdir}, salto", flush=True)
        return {'wmo': wmo, 'skipped': True, 'message': 'complete'}

    float_trend, df_hist = compute_float_trend(wmo, Hist_filtered_Profilelist, Dataset)
    save_df_report(out_meta, build_trend_report_row(float_trend), "Floats_trend_by_wmo.csv", tmp_report_dir=tmp_report_dir, wmo=wmo)

    if args_dict['make_plots'] and float_trend['DRIFT_CODE'] == 1:
        plot_drift_correction(wmo, df_hist, float_trend, out_meta)

    for p in Profilelist:
        outfile = get_outfile(p, outdir)
        doxy_exists = profile_doxy_exists.get(p.ID(), False)
        writing_mode = 'a' if profile_file_valid.get(p.ID(), False) else 'w'

        condition_to_write = not doxy_exists or args_dict['force']
        if not condition_to_write:
            continue

        metadata = Metadata(p._my_float.filename)
        doxy_algorithm(p, wmo, float_trend, df_hist, Dataset, outfile, metadata, writing_mode, out_meta, tmp_report_dir=tmp_report_dir)

    return {'wmo': wmo, 'skipped': False, 'message': 'processed'}


if __name__ == '__main__':
    input_file = args.update_file
    if input_file == 'NO_file':
        TI = TimeInterval(args.datestart, args.dateend, '%Y%m%d')
        R = Rectangle(-6,36,30,46)

        PROFILES_COR = bio_float.FloatSelector('DOXY', TI, R)
        wmo_list = bio_float.get_wmo_list(PROFILES_COR)
        wmo_list.sort()

        tmp_report_dir = Path(args.tmp_report_dir) if args.tmp_report_dir else OUT_META / 'tmp_reports'
        tmp_report_dir.mkdir(parents=True, exist_ok=True)

        args_dict = vars(args)
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(process_wmo, wmo, str(OUTDIR), str(OUT_META), args_dict): wmo for wmo in wmo_list}
            for future in as_completed(futures):
                wmo = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    print(f"{wmo} generated an exception: {exc}", flush=True)
                else:
                    print(f"{wmo} done: {result['message']}", flush=True)

        combine_temp_reports(tmp_report_dir, OUT_META, 'Floats_trend_by_wmo.csv')
        combine_temp_reports(tmp_report_dir, OUT_META, 'Floats_accepted.csv')
        combine_temp_reports(tmp_report_dir, OUT_META, 'Floats_rejected.csv')
        combine_temp_reports(tmp_report_dir, OUT_META, 'DataMode_and_Saturation_rejection_doxy.csv')

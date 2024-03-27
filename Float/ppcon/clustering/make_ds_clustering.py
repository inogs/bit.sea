import os
import torch
import numpy as np
import pandas as pd
import netCDF4 as nc
from discretization import dict_max_pressure, dict_interval
from make_ds.preprocessing import *
import sys
#sys.path.append("/g100_scratch/userexternal/camadio0/PPCON/bit.sea/")
from commons.time_interval import TimeInterval

max_pres_nitrate = dict_max_pressure["NITRATE"]
interval_nitrate = dict_interval["NITRATE"]

max_pres_chla = dict_max_pressure["CHLA"]
interval_chla = dict_interval["CHLA"]

max_pres_BBP700 = dict_max_pressure["BBP700"]
interval_BBP700 = dict_interval["BBP700"]


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def discretize(pres, var, max_pres, interval):
    """
    function that take as input a profile, the corresponding pres and create a tensor
    the interval input represents the discretization scale
    the max_pres input represents the maximum pressure considered
    """
    size = int(max_pres / interval)
    discretization_pres = np.arange(0, max_pres, interval)

    out = torch.zeros(size)

    for i in range(size):
        pressure_discretize = discretization_pres[i]
        idx = find_nearest(pres, pressure_discretize)
        if idx is None:
            return None
        out[i] = torch.from_numpy(np.asarray(var[idx]))

    return out


def read_date_time(date_time):
    year = int(date_time[0:4])
    month = int(date_time[4:6])
    month = month - 1
    day = int(date_time[6:8])

    day_total = month * 30 + day
    day_rad = day_total * 2 * np.pi / 365

    return year, day_rad


def make_dict_single_float(path, date_time):
    if not os.path.exists(path):
        return None

    year, day_rad = read_date_time(date_time)

    try:
        ds = nc.Dataset(path)  # Select sample
    except Exception as error:
        print('Caught this error: ' + repr(error))
        return dict()

    lat = ds["LATITUDE"][:].data[0]
    lon = ds["LONGITUDE"][:].data[0]
    # pres_df = ds["PRES"][:].data[:]
    psal_df = ds["PSAL"][:].data[:]
    pres_psal_df = ds["PRES_PSAL"][:].data[:]

    temp_df = ds["TEMP"][:].data[:]
    pres_temp_df = ds["PRES_TEMP"][:].data[:]

    if "DOXY" not in ds.variables.keys():
        return dict()
    doxy_df = ds["DOXY"][:].data[:]
    pres_doxy_df = ds["PRES_DOXY"][:].data[:]

    temp = discretize(pres_temp_df, temp_df, max_pres_nitrate, interval_nitrate)
    psal = discretize(pres_psal_df, psal_df, max_pres_nitrate, interval_nitrate)
    doxy = discretize(pres_doxy_df, doxy_df, max_pres_nitrate, interval_nitrate)

    if "NITRATE" not in ds.variables.keys():
        nitrate = torch.zeros(int(max_pres_nitrate / interval_nitrate))
    else:
        nitrate_df = ds["NITRATE"][:].data[:]
        pres_nitrate_df = ds["PRES_NITRATE"][:].data[:]
        nitrate = discretize(pres_nitrate_df, nitrate_df, max_pres_nitrate, interval_nitrate)

    if "CHLA" not in ds.variables.keys():
        chla = torch.zeros(int(max_pres_chla / interval_chla))
    else:
        chla_df = ds["CHLA"][:].data[:]
        pres_chla_df = ds["PRES_CHLA"][:].data[:]
        if missing_extreme_test(pres_variable_df=pres_chla_df, min_extreme=10,
                                max_extreme=180) and counting_measurement_test(variable_df=chla_df, tradeoff=50):
            chla = discretize(pres_chla_df, chla_df, max_pres_chla, interval_chla)
        else:
            return dict()

    if "BBP700" not in ds.variables.keys():
        BBP700 = torch.zeros(int(max_pres_BBP700 / interval_BBP700))
    else:
        BBP700_df = ds["BBP700"][:].data[:]
        pres_BBP700_df = ds["PRES_BBP700"][:].data[:]
        if missing_extreme_test(pres_variable_df=pres_BBP700_df, min_extreme=30,
                                max_extreme=180) and counting_measurement_test(variable_df=BBP700_df, tradeoff=50):
            BBP700_df = BBP700_df * 1000
            BBP700 = discretize(pres_BBP700_df, BBP700_df, max_pres_BBP700, interval_BBP700)
        else:
            return dict()

    name_float = path[8:-3]
    if temp is None or psal is None or doxy is None or nitrate is None or chla is None or BBP700 is None:
        return dict()
    dict_float = {name_float: [year, day_rad, lat, lon, temp, psal, doxy, nitrate, chla, BBP700]}

    return dict_float

"""
amadio new 
"""
def col_to_dt(df,name_col_object, name_col_datetime64):
    """ from object to datetime --> name_col
    """
    df[name_col_datetime64] = pd.to_datetime(df[name_col_object]).dt.date
    df[name_col_datetime64] = pd.to_datetime(df[name_col_datetime64])
    return(df)
"""
amadio modified inpit dir
"""
def make_dataset(path_float_index, INDIR):
    """ path_float_index  =  DIFF_floats_20240202.txt or Float_index.txt
        INDIR             = ONLINE_REPO
    """ 
 
    if path_float_index.split("/")[-1].startswith('DIFF_floats'):
       print(path_float_index)
       print(path_float_index)
       print(path_float_index)
       print(path_float_index)
       col           = ['filename','time','lat','lon','type','lenght','if','vars','vartype','updated']
       df            = pd.read_csv(path_float_index,header=None)
       df.columns=col
       # adjust the file diff as a floatindex.txt from the name to the datetype
       df.filename =  df.filename.str.replace("coriolis/", "") # remove coriolis
       df.filename =  df.filename.str.replace("profiles/", "") # remove profile
       df.time     =  df['time'].apply(str)
       df.time     =  df.time.str[:8]   # ony date as in float index 
       name_list     = df.filename.tolist()
       datetime_list = df.time.tolist()
       
    else: #update_file.split("/")[-1].startswith('Float_Index'):
       name_list     = pd.read_csv(path_float_index, header=None).to_numpy()[:, 0].tolist()
       datetime_list = pd.read_csv(path_float_index, header=None).to_numpy()[:, 3].tolist() 
    dict_ds_accepted = dict()

    for i in range(np.size(name_list)):
        path =INDIR+ "/SUPERFLOAT/" + name_list[i]
        if not os.path.exists(path):
            continue
        date_time = datetime_list[i]
        dict_single_float = make_dict_single_float(path, date_time)

        dict_ds_accepted = {**dict_ds_accepted, **dict_single_float}
    
    return dict_ds_accepted  # , dict_ds_removed


def make_pandas_df(path_float_index, path_clust_csv ,INDIR):
    """ path_float_index  =  DIFF_floats_20240202.txt or Float_index.txt
        path_clust_csv    =  ONLINE_REPO/clustering/ds_sf_clustering.csv
        INDIR             =  ONLINE_REPO
    """
    dict_ds_accepted = make_dataset(path_float_index, INDIR)
    pd_ds_accepted = pd.DataFrame(dict_ds_accepted,
                                  index=['year', 'day_rad', 'lat', 'lon', 'temp', 'psal', 'doxy', 'nitrate', 'chla',
                                         'BBP700'])
    pd_ds_accepted.to_csv( path_clust_csv  )
    return


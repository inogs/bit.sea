import os

import torch
import numpy as np
import pandas as pd
import netCDF4 as nc


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def find_nearest_attention(array, value, var):
    fill_values = [99999.0000]
    check_value = 99999.0000
    check_loop = 0
    excptIndx = []
    array = np.asarray(array)
    while check_value in fill_values and check_loop < len(array):
        differences = np.abs(array - value)

        m = np.zeros(differences.size, dtype=bool)
        m[excptIndx] = True
        differences_mask = np.ma.array(differences, mask=m)

        idx = np.argmin(differences_mask)

        check_value = var[idx]
        check_loop += 1
        excptIndx.append(idx)
    if check_loop == len(array):
        return None
    return idx


def discretize(pres, var, max_pres=2000, interval=10):
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
        idx = find_nearest_attention(pres, pressure_discretize, var)
        if idx is None:
            return None
        out[i] = torch.from_numpy(np.asarray(var[idx]))

    return out


def read_date_time(date_time):
    year = int(date_time[0:4])
    month = int(date_time[4:6])
    month = month - 1
    day = int(date_time[6:8])

    day_total = month + day
    day_rad = day_total * 2 * np.pi / 365

    return year, day_rad


def make_dict_single_float(path, date_time):
    if not os.path.exists(path):
        return None

    year, day_rad = read_date_time(date_time)

    ds = nc.Dataset(path)
    lat = ds["LATITUDE"][:].data[0]
    lon = ds["LONGITUDE"][:].data[0]
    pres_df = ds["PRES"][:].data[0]
    psal_df = ds["PSAL"][:].data[0]
    temp_df = ds["TEMP"][:].data[0]
    if "DOXY" not in ds.variables.keys():
        return dict()
    doxy_df = ds["DOXY"][:].data[0]
    if "NITRATE" not in ds.variables.keys():
        return dict()
    nitrate_df = ds["NITRATE"][:].data[0]

    temp = discretize(pres_df, temp_df)
    psal = discretize(pres_df, psal_df)
    doxy = discretize(pres_df, doxy_df)
    nitrate = discretize(pres_df, nitrate_df)
    name_float = path[8:-3]
    if temp is None or psal is None or doxy is None or nitrate is None:
        return dict()
    dict_float = {name_float: [year, day_rad, lat, lon, temp, psal, doxy, nitrate]}
    return dict_float


def make_dataset(path_float_index):
    name_list = pd.read_csv(path_float_index, header=None).to_numpy()[:, 0].tolist()
    datetime_list = pd.read_csv(path_float_index, header=None).to_numpy()[:, 3].tolist()

    dict_ds = dict()

    for i in range(np.size(name_list)):
        path = "CORIOLIS/" + name_list[i]
        if not os.path.exists(path):
            continue
        date_time = datetime_list[i]
        dict_single_float = make_dict_single_float(path, date_time)
        dict_ds = {**dict_ds, **dict_single_float}

    return dict_ds


def make_pandas_df(path_float_index):
    dict_ds = make_dataset(path_float_index)
    pd_ds = pd.DataFrame(dict_ds, index=['year', 'day_rad', 'lat', 'lon', 'temp', 'psal', 'doxy', 'nitrate'])

    pd_ds.to_csv(os.getcwd() + '/float_ds.csv')
    return


def make_toy_dataset(path_float_index):
    name_list = pd.read_csv(path_float_index, header=None).to_numpy()[:, 0].tolist()
    datetime_list = pd.read_csv(path_float_index, header=None).to_numpy()[:, 3].tolist()

    dict_ds = dict()

    for i in range(int(np.size(name_list)/20)):
        path = "FLOAT_BIO/data/" + name_list[i]
        if not os.path.exists(path):
            continue
        date_time = datetime_list[i]
        dict_single_float = make_dict_single_float(path, date_time)
        dict_ds = {**dict_ds, **dict_single_float}

    return dict_ds


def make_pandas_toy_df(path_float_index):
    dict_ds = make_toy_dataset(path_float_index)
    pd_ds = pd.DataFrame(dict_ds, index=['year', 'day_rad', 'lat', 'lon', 'temp', 'psal', 'doxy', 'nitrate'])

    pd_ds.to_csv(os.getcwd() + '/toy_ds.csv')
    return

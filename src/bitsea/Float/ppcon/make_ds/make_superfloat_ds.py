import os
import random

import torch
import numpy as np
import pandas as pd
import netCDF4 as nc

from bitsea.Float.ppcon.utils import shuffle_dict
from bitsea.Float.ppcon.discretization import dict_max_pressure, dict_interval
from bitsea.Float.ppcon.make_ds.preprocessing import *


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


def make_dict_single_float(path, date_time, variable):
    if not os.path.exists(path):
        return None

    year, day_rad = read_date_time(date_time)

    if year >= 2021:
        return dict(), 0

    try:
        ds = nc.Dataset(path)  # Select sample
    except Exception as error:
        print('Caught this error: ' + repr(error))
        return dict(), 0

    lat = ds["LATITUDE"][:].data[0]
    lon = ds["LONGITUDE"][:].data[0]
    # pres_df = ds["PRES"][:].data[:]
    psal_df = ds["PSAL"][:].data[:]
    pres_psal_df = ds["PRES_PSAL"][:].data[:]

    temp_df = ds["TEMP"][:].data[:]
    pres_temp_df = ds["PRES_TEMP"][:].data[:]

    if "DOXY" not in ds.variables.keys():
        return dict(), 0
    doxy_df = ds["DOXY"][:].data[:]
    pres_doxy_df = ds["PRES_DOXY"][:].data[:]

    if variable not in ds.variables.keys():
        return dict(), 0
    variable_df = ds[f"{variable}"][:].data[:]
    pres_variable_df = ds[f"PRES_{variable}"][:].data[:]

    if variable == "NITRATE":
        flag_missing_extreme = 1
        flag_counting_measurement = 1

    if variable == "CHLA":
        flag_missing_extreme = missing_extreme_test(pres_variable_df=pres_variable_df, min_extreme=10, max_extreme=180)
        flag_counting_measurement = counting_measurement_test(variable_df=variable_df, tradeoff=50)

    if variable == "BBP700":
        flag_missing_extreme = missing_extreme_test(pres_variable_df=pres_variable_df, min_extreme=30, max_extreme=180)
        flag_counting_measurement = counting_measurement_test(variable_df=variable_df, tradeoff=50)
        variable_df = variable_df * 1000

    max_pres = dict_max_pressure[variable]
    interval = dict_interval[variable]

    temp = discretize(pres_temp_df, temp_df, max_pres, interval)
    psal = discretize(pres_psal_df, psal_df, max_pres, interval)
    doxy = discretize(pres_doxy_df, doxy_df, max_pres, interval)

    variable = discretize(pres_variable_df, variable_df, max_pres, interval)
    name_float = path[8:-3]
    if temp is None or psal is None or doxy is None or variable is None:
        return dict(), 0

    dict_float = {name_float: [year, day_rad, lat, lon, temp, psal, doxy, variable]}
    flag = flag_missing_extreme * flag_counting_measurement  # if at least one of the flag is 0 then the final flag
    # is zero

    return dict_float, flag


def make_dataset(path_float_index, variable):
    name_list = pd.read_csv(path_float_index, header=None).to_numpy()[:, 0].tolist()
    datetime_list = pd.read_csv(path_float_index, header=None).to_numpy()[:, 3].tolist()

    dict_ds_accepted = dict()
    dict_ds_removed = dict()

    for i in range(np.size(name_list)):
        path = os.getcwd() + "/ds/SUPERFLOAT/" + name_list[i]
        if not os.path.exists(path):
            continue
        date_time = datetime_list[i]
        dict_single_float, flag = make_dict_single_float(path, date_time, variable)
        # dict_single_float = make_dict_single_float(path, date_time, variable)

        # dict_ds_accepted = {**dict_ds_accepted, **dict_single_float}
        if flag == 1:
            dict_ds_accepted = {**dict_ds_accepted, **dict_single_float}
        if flag == 0:
            dict_ds_removed = {**dict_ds_removed, **dict_single_float}

    return dict_ds_accepted, dict_ds_removed


def make_pandas_df(path_float_index, variable):
    dict_ds_accepted, dict_ds_removed = make_dataset(path_float_index, variable)
    dict_ds_accepted = shuffle_dict(dict_ds_accepted)
    # dict_ds_accepted = make_dataset(path_float_index, variable)

    pd_ds_removed = pd.DataFrame(dict_ds_removed,
                                 index=['year', 'day_rad', 'lat', 'lon', 'temp', 'psal', 'doxy', variable])

    train_frac = 0.8
    train_size = int(train_frac * len(dict_ds_accepted))
    val_size = len(dict_ds_accepted) - train_size

    dict_ds_accepted_train = dict_ds_accepted.copy()
    dict_ds_accepted_test = dict(list(dict_ds_accepted.items())[:val_size])
    for key_test in dict_ds_accepted_test.keys():
        dict_ds_accepted_train.pop(key_test, None)

    pd_ds_accepted_train = pd.DataFrame(dict_ds_accepted_train,
                                        index=['year', 'day_rad', 'lat', 'lon', 'temp', 'psal', 'doxy', variable])
    pd_ds_accepted_test = pd.DataFrame(dict_ds_accepted_test,
                                       index=['year', 'day_rad', 'lat', 'lon', 'temp', 'psal', 'doxy', variable])

    # train_ds_accepted, test_ds_accepted = train_test_split(pd_ds_accepted, test_size=0.2, shuffle=True)

    pd_ds_accepted_train.to_csv(os.getcwd() + f'/ds/{variable}/float_ds_sf_train.csv')
    pd_ds_accepted_test.to_csv(os.getcwd() + f'/ds/{variable}/float_ds_sf_test.csv')

    pd_ds_removed.to_csv(os.getcwd() + f'/ds/{variable}/float_ds_sf_removed.csv')

    return


def make_toy_dataset(path_float_index, variable):
    name_list = pd.read_csv(path_float_index, header=None).to_numpy()[:, 0].tolist()
    datetime_list = pd.read_csv(path_float_index, header=None).to_numpy()[:, 3].tolist()

    dict_ds = dict()

    for i in range(int(np.size(name_list) / 20)):
        path = os.getcwd() + "/ds/SUPERFLOAT/" + name_list[i]
        if not os.path.exists(path):
            continue
        date_time = datetime_list[i]
        dict_single_float = make_dict_single_float(path, date_time, variable)
        dict_ds = {**dict_ds, **dict_single_float}

    return dict_ds


def make_pandas_toy_df(path_float_index, variable):
    dict_ds = make_toy_dataset(path_float_index, variable)
    pd_ds = pd.DataFrame(dict_ds, index=['year', 'day_rad', 'lat', 'lon', 'temp', 'psal', 'doxy', variable])

    pd_ds.to_csv(os.getcwd() + f'/ds/{variable}/toy_ds_sf.csv')
    return

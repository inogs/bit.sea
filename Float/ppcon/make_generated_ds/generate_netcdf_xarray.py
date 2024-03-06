import os

import xarray as xr
import numpy as np
import pandas as pd
import netCDF4 as nc
import torch
from torch.utils.data import DataLoader

from discretization import *
from utils import upload_and_evaluate_model, get_output
from dataset_with_float_names import FloatDataset


dict_models = {
    "NITRATE": ["2023-04-04_", 50],
    "CHLA": ["2023-03-29", 150],
    "BBP700": ["2023-03-29", 100]
}

date_nitrate = dict_models["NITRATE"][0]
date_chla = dict_models["CHLA"][0]
date_BBP700 = dict_models["BBP700"][0]

dir_model_nitrate = os.getcwd() + f"/../results/NITRATE/{date_nitrate}/model"
dir_model_chla = os.getcwd() + f"/../results/CHLA/{date_chla}/model"
dir_model_BBP700 = os.getcwd() + f"/../results/BBP700/{date_BBP700}/model"

# Upload the input information
info_nitrate = pd.read_csv(os.getcwd() + f"/../results/NITRATE/{date_nitrate}/info.csv")
info_chla = pd.read_csv(os.getcwd() + f"/../results/CHLA/{date_chla}/info.csv")
info_BBP700 = pd.read_csv(os.getcwd() + f"/../results/BBP700/{date_BBP700}/info.csv")

# Upload and evaluate the model
model_day_nitrate, model_year_nitrate, model_lat_nitrate, model_lon_nitrate, model_nitrate = upload_and_evaluate_model(
    dir_model=dir_model_nitrate, info_model=info_nitrate, ep=dict_models["NITRATE"][1])
model_day_chla, model_year_chla, model_lat_chla, model_lon_chla, model_chla = upload_and_evaluate_model(
    dir_model=dir_model_chla, info_model=info_chla, ep=dict_models["CHLA"][1])
model_day_BBP700, model_year_BBP700, model_lat_BBP700, model_lon_BBP700, model_BBP700 = upload_and_evaluate_model(
    dir_model=dir_model_BBP700, info_model=info_BBP700, ep=dict_models["BBP700"][1])

# create the new enhanced directory
path_superfloat_ppcon = f"/home/gpietropolli/Desktop/canyon-float/ds/SUPERFLOAT_PPCon/"
if not os.path.exists(path_superfloat_ppcon):
    os.mkdir(path_superfloat_ppcon)

# read the float number + sampling number from the dataset
path_df = f"/home/gpietropolli/Desktop/canyon-float/ds/clustering/ds_sf_clustering.csv"
dataset = FloatDataset(path_df)
my_ds = DataLoader(dataset, shuffle=True)

# information related to the pres of PPCon generated variables
pres_nitrate = np.arange(0, dict_max_pressure["NITRATE"], dict_interval["NITRATE"])
pres_chla = np.arange(0, dict_max_pressure["CHLA"], dict_interval["CHLA"])
pres_BBP700 = np.arange(0, dict_max_pressure["BBP700"], dict_interval["BBP700"])

for sample in my_ds:
    year, day_rad, lat, lon, temp, psal, doxy, nitrate, chla, BBP700, name_float = sample
    sample_prediction = sample[:-1]

    name_file = name_float[0]
    main_dir = name_float[0][2:-4]
    print(f"generating profiles {name_file}")
    # create the directory in the enhanced dataset direcotry
    path_saving_nc = path_superfloat_ppcon + f"{main_dir}/"
    if not os.path.exists(path_saving_nc):
        os.mkdir(path_saving_nc)

    # open the netcfd file in the "ds/SUPERFLOAT" directory
    path_superfloat = f"/home/gpietropolli/Desktop/canyon-float/ds/SUPERFLOAT/{main_dir}/{name_file}.nc"
    ds = xr.open_dataset(path_superfloat)  # ds = nc.Dataset(path_superfloat)

    # create the dimension relative to the PPCon prediction
    # PPCon_dim = ds.createDimension('PPCon_dim', 200)

    # apply CNN model and save the updated netcdf file
    if torch.count_nonzero(nitrate) == 0:
        nitrate = get_output(sample=sample_prediction, model_day=model_day_nitrate, model_year=model_year_nitrate,
                             model_lat=model_lat_nitrate, model_lon=model_lon_nitrate, model=model_nitrate)
        nitrate = nitrate.detach()

    nitrate = torch.squeeze(nitrate)
    nitrate = nitrate.numpy()

    ds["NITRATE_PPCon"] = xr.DataArray(nitrate)
    ds["NITRATE_PPCon_PRES"] = xr.DataArray(pres_nitrate)

    if torch.count_nonzero(chla) == 0:
        chla = get_output(sample=sample_prediction, model_day=model_day_chla, model_year=model_year_chla,
                          model_lat=model_lat_chla, model_lon=model_lon_chla, model=model_chla)
        chla = chla.detach()

    chla = torch.squeeze(chla)
    chla = chla.numpy()

    ds["CHLA_PPCon"] = xr.DataArray(chla)
    ds["CHLA_PPCon_PRES"] = xr.DataArray(pres_chla)

    if torch.count_nonzero(BBP700) == 0:
        BBP700 = get_output(sample=sample_prediction, model_day=model_day_BBP700, model_year=model_year_BBP700,
                            model_lat=model_lat_BBP700, model_lon=model_lon_BBP700, model=model_BBP700)
        BBP700 = BBP700.detach()
        BBP700 = BBP700 / 1000

    BBP700 = torch.squeeze(BBP700)
    BBP700 = BBP700.numpy()

    ds["BBP700_PPCon"] = xr.DataArray(BBP700)
    ds["BBP700_PPCon_PRES"] = xr.DataArray(pres_BBP700)

    # save the new generated nc file in the enhanced directory
    ds.to_netcdf(f"{path_saving_nc}/{name_file}.nc")







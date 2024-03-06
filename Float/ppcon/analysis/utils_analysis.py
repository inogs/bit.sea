import os
import random

import numpy as np
import pandas as pd
import torch

from torch.nn.functional import mse_loss
from torch.utils.data import DataLoader

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
#from mpl_toolkits.basemap import Basemap

from discretization import dict_max_pressure
from dataset import FloatDataset
from utils import upload_and_evaluate_model, from_day_rad_to_day
from other_methods.suazade import get_suazade_profile
from other_methods.gloria import get_gloria_profile


sns.set_theme(context='paper', style='whitegrid', font='sans-serif', font_scale=1.5,
              color_codes=True, rc=None)

dict_unit_measure = {"NITRATE": r"$mmol \ m^{-3}$",
                     "CHLA": r"$mg \ m^{-3}$",
                     "BBP700": r"$m^{-1}$"}

dict_var_name = {"NITRATE": "Nitrate",
                 "CHLA": "Chlorophyll",
                 "BBP700": "BBP700"}

pal = sns.color_palette("magma")

dict_color = {'NWM': pal[0], 'SWM': pal[1], 'TYR': pal[3], 'ION': pal[4], 'LEV': pal[5]}


def count_samples(variable):
    path_ds = os.getcwd() + f"/../ds/{variable}/"
    ds_train = FloatDataset(path_ds + "float_ds_sf_train.csv")
    ds_test = FloatDataset(path_ds + "float_ds_sf_test.csv")
    ds_removed = FloatDataset(path_ds + "float_ds_sf_removed.csv")
    return len(ds_train), len(ds_test), len(ds_removed)


def get_reconstruction(variable, date_model, epoch_model, mode):
    # Upload the input ds
    path_float = os.getcwd() + f"/../ds/{variable}/float_ds_sf_{mode}.csv"
    if mode == "all":
        path_float = os.getcwd() + f"/../ds/{variable}/float_ds_sf.csv"
    dataset = FloatDataset(path_float)
    ds = DataLoader(dataset, shuffle=True)

    dir_model = os.getcwd() + f"/../results/{variable}/{date_model}/model"
    info = pd.read_csv(os.getcwd() + f"/../results/{variable}/{date_model}/info.csv")

    # Upload and evaluate the model
    model_day, model_year, model_lat, model_lon, model = upload_and_evaluate_model(dir_model=dir_model, info_model=info,
                                                                                   ep=epoch_model)

    lat_list = list()
    lon_list = list()
    day_rad_list = list()
    measured_var_list = list()
    generated_var_list = list()

    for sample in ds:
        year, day_rad, lat, lon, temp, psal, doxy, measured_var = sample

        output_day = model_day(day_rad.unsqueeze(1))
        output_year = model_year(year.unsqueeze(1))
        output_lat = model_lat(lat.unsqueeze(1))
        output_lon = model_lon(lon.unsqueeze(1))

        output_day = torch.transpose(output_day.unsqueeze(0), 0, 1)
        output_year = torch.transpose(output_year.unsqueeze(0), 0, 1)
        output_lat = torch.transpose(output_lat.unsqueeze(0), 0, 1)
        output_lon = torch.transpose(output_lon.unsqueeze(0), 0, 1)
        temp = torch.transpose(temp.unsqueeze(0), 0, 1)
        psal = torch.transpose(psal.unsqueeze(0), 0, 1)
        doxy = torch.transpose(doxy.unsqueeze(0), 0, 1)

        x = torch.cat((output_day, output_year, output_lat, output_lon, temp, psal, doxy), 1)
        generated_var = model(x.float())
        generated_var = generated_var.detach()  # torch.squeeze????

        lat_list.append(lat.item())
        lon_list.append(lon.item())
        day_rad_list.append(day_rad.item())

        if variable == "NITRATE":
            generated_var = torch.squeeze(generated_var)[:-10]
            measured_var = torch.squeeze(measured_var)[:-10]
        if variable == "BBP700":
            generated_var = torch.squeeze(generated_var) / 1000
            measured_var = torch.squeeze(measured_var) / 1000
        else:
            generated_var = torch.squeeze(generated_var)
            measured_var = torch.squeeze(measured_var)

        generated_var_list.append(generated_var)
        measured_var_list.append(measured_var)

    return lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list


def compute_efficency(generated_var_list, measured_var_list, variance):
    list_mse = []
    loss_mse = torch.nn.MSELoss(reduction='none')

    for k in range(len(generated_var_list)):
        list_mse.append(loss_mse(generated_var_list[k], measured_var_list[k]))
    dim_profile = int(generated_var_list[0].size(dim=0))
    tensor_sum_mse = torch.zeros(dim_profile)
    element = 0
    for el in range(len(list_mse)):
        tensor_sum_mse += list_mse[el]
        element += 1
    tensor_mse = tensor_sum_mse / element
    tensor_rmse = torch.sqrt(tensor_mse)

    efficency = 1 - (tensor_rmse / variance)

    return efficency

import os

import numpy as np
import pandas as pd
import torch

import matplotlib.pyplot as plt

# from conv1med import Conv1dMed
from architecture.conv1med_dp import Conv1dMed
from architecture.mlp import MLPDay, MLPYear, MLPLon, MLPLat

from discretization import dict_max_pressure, dict_interval
import scipy.ndimage


def plot_profiles(ds, dir, variable, ep=100):

    dir_model = dir + "/model"
    dir_profile = dir + "/profile"
    if not os.path.exists(dir_profile):
        os.mkdir(path=dir_profile)

    dir_info = dir + "/info.csv"
    info_df = pd.read_csv(dir_info)
    dp_rate = info_df['dp_rate'].item()

    max_pres = dict_max_pressure[variable]
    interval = dict_interval[variable]

    # Path of the saved models
    path_model_day = dir_model + "/model_day_" + str(ep) + ".pt"
    path_model_year = dir_model + "/model_year_" + str(ep) + ".pt"
    path_model_lon = dir_model + "/model_lon_" + str(ep) + ".pt"
    path_model_lat = dir_model + "/model_lat_" + str(ep) + ".pt"
    path_model_conv = dir_model + "/model_conv_" + str(ep) + ".pt"

    # Upload and evaluate all the models necessary
    # """
    model_day = MLPDay()
    model_day.load_state_dict(torch.load(path_model_day))
    model_day.eval()

    model_year = MLPYear()
    model_year.load_state_dict(torch.load(path_model_year))
    model_year.eval()

    model_lat = MLPLat()
    model_lat.load_state_dict(torch.load(path_model_lat))
    model_lat.eval()

    model_lon = MLPLon()
    model_lon.load_state_dict(torch.load(path_model_lon))
    model_lon.eval()

    model = Conv1dMed(dp_rate=dp_rate)
    model.load_state_dict(torch.load(path_model_conv))
    model.eval()

    for year, day_rad, lat, lon, temp, psal, doxy, output_variable in ds:

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
        output_variable = torch.transpose(output_variable.unsqueeze(0), 0, 1)

        x = torch.cat((output_day, output_year, output_lat, output_lon, temp, psal, doxy), 1)
        output_test = model(x.float())

        depth_output = np.linspace(0, max_pres, len(output_test[0, 0, :].detach().numpy()))
        depth_variable = np.linspace(0, max_pres, len(output_variable[0, 0, :].detach().numpy()))

        if variable == "BBP700":
            output_variable = output_variable / 1000
            output_test = output_test / 1000

        plt.plot(output_test[0, 0, :].detach().numpy(), depth_output, label=f"generated {variable}")
        plt.plot(output_variable[0, 0, :].detach().numpy(), depth_variable, label=f"measured {variable}")
        plt.gca().invert_yaxis()

        plt.legend()
        plt.savefig(dir_profile + f"/profile_{year.item()}_{day_rad.item()}_{round(lat.item(), 2)}_{round(lon.item(), 2)}.png")
        # plt.show()
        plt.close()




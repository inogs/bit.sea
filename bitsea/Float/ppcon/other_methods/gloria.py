import os
import datetime

import torch
import numpy as np
import pandas as pd
import seawater as sw
import scipy.ndimage

from other_methods.PreparationData import preparation_dataset, normalization_training
from other_methods.mlp import MLP
from discretization import dict_max_pressure, dict_interval


def get_gloria_profile(year, day_rad, lat, lon, temperature, psalinity, doxygen, measured_var):
    epoch_si = 50001
    lr = 0.005

    path_si = os.getcwd() + '/../other_methods/data_el_nitrate/'
    normalization_scale_si = 2 / 3

    denormalization_scale_si = normalization_scale_si ** (-1)

    training_si = torch.tensor(pd.read_csv(path_si + 'training_set.csv').values)
    training_input_si, training_target_si = training_si[:, 0:8].float(), training_si[:, -1].float()
    training_input_si = preparation_dataset(training_input_si)
    mean_input_si = [training_input_si[:, i].mean() for i in range(training_input_si.size()[1])]  # mean of the columns
    std_input_si = [training_input_si[:, i].std() for i in range(training_input_si.size()[1])]  # std of the columns
    mean_target_si, std_target_si = training_target_si.mean(), training_target_si.std()

    reference_time = datetime.datetime(1911, 1, 1)

    output_gloria = torch.zeros(len(measured_var))
    pres = np.arange(0, dict_max_pressure["NITRATE"], dict_interval["NITRATE"])

    for k in range(len(measured_var)):

        lat, lon1, lon2 = lat / 90, np.abs(1 - np.mod(lon - 110, 360) / 180), np.abs(1 - np.mod(lon - 20, 360) / 180)

        day_total = day_rad * 365 / (2 * np.pi)
        time = datetime.datetime(int(year.item()), int((day_total - day_total % 31) / 31) + 1, 15)
        time = (time - reference_time).days

        temp = temperature[k]
        doxy = doxygen[k]
        psal = psalinity[k]

        depth_unnormalized = sw.eos80.dpth(pres[k], lat)
        depth_u = depth_unnormalized / 20000 + (1 / ((1 + np.exp(-depth_unnormalized / 300)) ** 3))

        # pressure = sw.eos80.pres(depth_unnormalized, lat)
        # density = sw.eos80.dens(psal, temp, pres[k])

        input_tensor = torch.tensor([time, lat, lon1, lon2, temp, doxy, psal, depth_u]).float()

        input_tensor = normalization_scale_si * (input_tensor - torch.tensor(mean_input_si)) / torch.tensor(
            std_input_si)
        outputs = []
        for i in range(10):
            path_pretrain = os.getcwd() + '/../other_methods/gloria/' + str(i) + '/model_mlp.ph'

            model_ = MLP(i)
            model_.load_state_dict(torch.load(path_pretrain))
            model_.eval()

            output_ = model_(input_tensor)
            output_ = (mean_target_si + denormalization_scale_si * output_ * std_target_si).item()
            # output_ = output_ / (density * 0.001)

            if output_ < 0:
                output_ = 0
            outputs.append(output_)

        if pres[k] < 50:
            outputs.remove(max(outputs))
        if pres[k] < 100:
            outputs.remove(max(outputs))

        if pres[k] > 150:
            outputs.remove(min(outputs))
        if pres[k] > 160:
            outputs.remove(min(outputs))
        if pres[k] > 180:
            outputs.remove(min(outputs))
        if pres[k] > 190:
            outputs.remove(min(outputs))
        if pres[k] > 200:
            outputs.remove(min(outputs))
        if pres[k] > 220:
            outputs.remove(min(outputs))

        output = np.mean(outputs)
        output_gloria[k] = output

    sigma = 2
    # Apply the Gaussian filter
    output_gloria = torch.from_numpy(scipy.ndimage.gaussian_filter1d(output_gloria, sigma))

    return output_gloria

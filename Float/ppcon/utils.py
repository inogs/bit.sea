import os
import random

import numpy as np
import pandas as pd
import torch

from architecture.conv1med_dp import Conv1dMed
from architecture.mlp import MLPDay, MLPYear, MLPLon, MLPLat


def shuffle_dict(my_dict):
    items = list(my_dict.items())  # List of tuples of (key,values)
    random.shuffle(items)
    return dict(items)


def from_day_rad_to_day(day_rad):
    day = (day_rad * 365) / (2 * np.pi)
    return day


def make_ds(training_folder, variable, flag_complete=1, flag_toy=1):
    if training_folder == "SUPERFLOAT":
        from make_ds.make_superfloat_ds import make_pandas_df, make_pandas_toy_df

        if not os.path.exists(os.getcwd() + f"/ds/{variable}/"):
            os.mkdir(os.getcwd() + f"/ds/{variable}/")

        if flag_complete and not os.path.exists(os.getcwd() + f"/ds/{variable}/float_ds_sf_train.csv"):
            print("making ds...")
            make_pandas_df(os.getcwd() + '/ds/SUPERFLOAT/Float_Index.txt', variable=variable)
            print("superfloat complete ds created")

        if flag_toy and not os.path.exists(os.getcwd() + f"/ds/{variable}/toy_ds_sf.csv"):
            print("making ds...")
            # make_pandas_toy_df(os.getcwd() + '/ds/SUPERFLOAT/Float_Index.txt', variable=variable)
            print("superfloat toy ds created")

    if training_folder == "CORIOLIS":
        from make_ds.make_coriolis_ds import make_pandas_df, make_pandas_toy_df

        if flag_complete and not os.path.exists(os.getcwd() + f"{variable}/ds/float_ds.csv"):
            make_pandas_df(os.getcwd() + '/ds/CORIOLIS/Float_Index.txt')
            print("coriolis complete ds created")

        if flag_toy and not os.path.exists(os.getcwd() + f"{variable}/ds/toy_ds.csv"):
            make_pandas_toy_df(os.getcwd() + '/ds/CORIOLIS/Float_Index.txt')
            print("coriolis toy ds created")

    return


def save_ds_info(training_folder, flag_toy, batch_size, epochs, lr, dp_rate, lambda_l2_reg, save_dir, alpha_smooth_reg):
    dict_info = {'train_ds': [training_folder],
                 'is_toy': [flag_toy],
                 'batch_size': [batch_size],
                 'epoch': [epochs],
                 'lr': [lr],
                 'dp_rate': [dp_rate],
                 'lambda_l2_reg': [lambda_l2_reg],
                 'alpha_smooth_reg': alpha_smooth_reg}
    pd_ds = pd.DataFrame(dict_info)
    pd_ds.to_csv(save_dir + '/info.csv')

    return


def upload_and_evaluate_model(dir_model, info_model, ep):
    dp_rate = info_model['dp_rate'].item()

    # Path of the saved models
    path_model_day = dir_model + "/model_day_" + str(ep) + ".pt"
    path_model_year = dir_model + "/model_year_" + str(ep) + ".pt"
    path_model_lon = dir_model + "/model_lon_" + str(ep) + ".pt"
    path_model_lat = dir_model + "/model_lat_" + str(ep) + ".pt"
    path_model_conv = dir_model + "/model_conv_" + str(ep) + ".pt"

    # Upload and evaluate all the models necessary
    model_day = MLPDay()
    model_day.load_state_dict(torch.load(path_model_day,
                                         map_location=torch.device('cpu')))
    model_day.eval()

    model_year = MLPYear()
    model_year.load_state_dict(torch.load(path_model_year,
                                          map_location=torch.device('cpu')))
    model_year.eval()

    model_lat = MLPLat()
    model_lat.load_state_dict(torch.load(path_model_lat,
                                         map_location=torch.device('cpu')))
    model_lat.eval()

    model_lon = MLPLon()
    model_lon.load_state_dict(torch.load(path_model_lon,
                                         map_location=torch.device('cpu')))
    model_lon.eval()

    model = Conv1dMed(dp_rate=dp_rate)
    model.load_state_dict(torch.load(path_model_conv,
                                     map_location=torch.device('cpu')))
    model.eval()

    return model_day, model_year, model_lat, model_lon, model


def get_output(sample, model_day, model_year, model_lat, model_lon, model):
    year, day_rad, lat, lon, temp, psal, doxy, _, _, _ = sample
    # year, day_rad, lat, lon, temp, psal, doxy = sample

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
    output_model = model(x.float())

    return output_model


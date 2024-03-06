import os

import torch

from dataset import FloatDataset


def get_reconstruction(variable, date_model, epoch_model, mode):
    # Upload the input ds
    path_float = f"/home/gpietropolli/Desktop/canyon-float/ds/{variable}/float_ds_sf_{mode}.csv"
    if mode == "all":
        path_float = f"/home/gpietropolli/Desktop/canyon-float/ds/{variable}/float_ds_sf.csv"
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
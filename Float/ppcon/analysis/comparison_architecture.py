import os
import scipy

from utils_analysis import *


def moving_average(data, window_size):
    if window_size % 2 == 0:
        window_size += 1  # Ensure window size is odd for symmetry
    pad_size = window_size // 2
    padded_data = np.pad(data, pad_size, mode='edge')
    cumsum_vec = np.cumsum(np.insert(padded_data, 0, 0))
    return (cumsum_vec[window_size:] - cumsum_vec[:-window_size]) / window_size


def reconstruction_profiles(variable, date_model, epoch_model, mode):
    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/fig/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)
    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/fig/profile_{mode}_{epoch_model}/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    number_samples = len(generated_var_list)

    for index_sample in range(number_samples):
        lat = lat_list[index_sample]
        lon = lon_list[index_sample]
        generated_profile = generated_var_list[index_sample]
        measured_profile = measured_var_list[index_sample]

        max_pres = dict_max_pressure[variable]
        depth = np.linspace(0, max_pres, len(generated_profile.detach().numpy()))
        plt.figure(figsize=(6, 7))

        measured_profile = moving_average(measured_profile.detach().numpy(), 3)
        generated_profile = moving_average(generated_profile.detach().numpy(), 3)
        # generated_profile = moving_average(generated_profile.detach().numpy(), 3)

        plt.plot(measured_profile, depth, lw=3, color="#2CA02C", label=f"Measured")
        plt.plot(generated_profile, depth, lw=3, linestyle=(0, (3, 1, 1, 1)), color="#1F77B4", label=f"PPCon")
        plt.gca().invert_yaxis()

        plt.xlabel(f"{dict_var_name[variable]} [{dict_unit_measure[variable]}]")
        plt.ylabel(r"Depth [$m$]")

        if variable == "BBP700":
            ax = plt.gca()
            x_labels = ax.get_xticks()
            ax.set_xticklabels(['{:,.0e}'.format(x) for x in x_labels])

        plt.legend()
        plt.tight_layout()

        plt.legend()
        plt.savefig(f"{path_analysis}profile_{round(lat, 2)}_{round(lon, 2)}.png") #, dpi=1200)
        plt.close()

    return


def reconstruction_profile_MLP(variable, date_model, epoch_model, mode):
    path_analysis = os.getcwd() + f"/../results/NITRATE/{date_model}/fig/comparison-def"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    path_float = os.getcwd() + f"/../ds/NITRATE/float_ds_sf_{mode}.csv"
    if mode == "all":
        path_float = os.getcwd() + f"/../ds/NITRATE/float_ds_sf.csv"
    dataset = FloatDataset(path_float)
    ds = DataLoader(dataset, shuffle=True)

    dir_model = os.getcwd() + f"/../results/NITRATE/{date_model}/model"
    info = pd.read_csv(os.getcwd() + f"/../results/NITRATE/{date_model}/info.csv")

    model_day, model_year, model_lat, model_lon, model = upload_and_evaluate_model(dir_model=dir_model, info_model=info,
                                                                                   ep=epoch_model)

    lat_list = list()
    lon_list = list()
    day_rad_list = list()
    measured_var_list = list()
    generated_var_list = list()
    gloria_var_list = list()
    canyon_var_list = list()

    sum_generated = torch.zeros(1, 1, 200)
    sum_measured = torch.zeros(1, 1, 200)
    sum_gloria = torch.zeros(1, 1, 200)
    sum_canyon = torch.zeros(1, 1, 200)

    loss_gloria = 0
    loss_canyon = 0
    loss_ppcon = 0

    number_profiles = 0

    for sample in ds:
        year, day_rad, lat, lon, temp, psal, doxy, measured_var = sample

        number_profiles += 1

        generated_gloria_var = get_gloria_profile(year, day_rad, lat, lon, torch.squeeze(temp), torch.squeeze(psal),
                                                  torch.squeeze(doxy), torch.squeeze(measured_var))
        gloria_var_list.append(generated_gloria_var)
        sum_gloria += generated_gloria_var
        loss_gloria += np.sqrt(mse_loss(generated_gloria_var, measured_var.squeeze()))
        print(f"MLP Pietropolli: {mse_loss(generated_gloria_var, measured_var.squeeze())}")

        generated_canyon_var = get_suazade_profile(year, day_rad, lat, lon, torch.squeeze(temp), torch.squeeze(psal),
                                                   torch.squeeze(doxy), torch.squeeze(measured_var))
        canyon_var_list.append(generated_canyon_var)
        sum_canyon += generated_canyon_var
        loss_canyon += np.sqrt(mse_loss(generated_canyon_var, measured_var.squeeze()))
        print(f"CANYON-Med: {mse_loss(generated_canyon_var, measured_var.squeeze())}")

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

        sigma = 2
        generated_var = torch.from_numpy(scipy.ndimage.gaussian_filter1d(generated_var, sigma))

        loss_ppcon += np.sqrt(mse_loss(generated_var.squeeze(), measured_var.squeeze()))
        print(f"PPCon after reg: {mse_loss(generated_var.squeeze(), measured_var.squeeze())}")
        print(f"MLP Pietropolli: {loss_gloria / number_profiles} \t CANYON-Med: {loss_canyon / number_profiles} \t PPCon: {loss_ppcon / number_profiles}")

        with open(os.getcwd() + '/../results.txt', 'a') as file:
            file.write(
                f"MLP Pietropolli: {loss_gloria / number_profiles} \t CANYON-Med: {loss_canyon / number_profiles} \t PPCon: {loss_ppcon / number_profiles}\n")

        lat_list.append(lat.item())
        lon_list.append(lon.item())
        day_rad_list.append(day_rad.item())

        generated_var_list.append(generated_var)
        sum_generated += generated_var
        measured_var_list.append(measured_var)
        sum_measured += measured_var

        depth = np.linspace(0, dict_max_pressure["NITRATE"], len(torch.squeeze(generated_var)))

        cut_sup = 3
        cut_inf = 10

        measured_var = torch.squeeze(measured_var)[cut_sup:-cut_inf]
        measured_var = moving_average(measured_var.numpy(), 5)

        plt.figure(figsize=(6, 7))

        plt.plot(measured_var, depth[cut_sup:-cut_inf], lw=3, color="#2CA02C",
                 label="Measured")
        plt.plot(torch.squeeze(generated_gloria_var)[cut_sup:-cut_inf], depth[cut_sup:-cut_inf], lw=3,
                 color="midnightblue", linestyle=(0, (5, 1)), label="MLP")
        plt.plot(torch.squeeze(generated_canyon_var)[cut_sup:-cut_inf], depth[cut_sup:-cut_inf], lw=3,
                 color="blue", linestyle=(0, (5, 1)), label="CANYON-Med")
        plt.plot(torch.squeeze(generated_var)[cut_sup:-cut_inf], depth[cut_sup:-cut_inf], lw=3,
                 color="#1F77B4", linestyle=(0, (3, 1, 1, 1)), label="PPCon")
        plt.gca().invert_yaxis()

        plt.xlabel(f"{variable} [{dict_unit_measure[variable]}]")
        plt.ylabel(r"depth [$m$]")

        plt.legend()
        plt.tight_layout()

        # plt.show()
        plt.savefig(f"{path_analysis}/method_comparison_{round(lat.item(), 2)}_{round(lon.item(), 2)}.png")  # , dpi=600)
        plt.close()

    return
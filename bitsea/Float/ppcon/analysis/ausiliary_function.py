from utils_analysis import *


def get_reconstruction_comparison(variable, season, date_model, epoch_model, mode):
    # Create savedir
    path_analysis = os.getcwd() + f"/../results/NITRATE/{date_model}/fig/comp"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)
    # Upload the input ds

    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}

    path_float = f"/home/gpietropolli/Desktop/canyon-float/ds/NITRATE/float_ds_sf_{mode}.csv"
    if mode == "all":
        path_float = f"/home/gpietropolli/Desktop/canyon-float/ds/NITRATE/float_ds_sf.csv"
    dataset = FloatDataset(path_float)
    ds = DataLoader(dataset, shuffle=True)

    dir_model = os.getcwd() + f"/../results/NITRATE/{date_model}/model"
    info = pd.read_csv(os.getcwd() + f"/../results/NITRATE/{date_model}/info.csv")

    # Upload and evaluate the model
    model_day, model_year, model_lat, model_lon, model = upload_and_evaluate_model(dir_model=dir_model, info_model=info,
                                                                                   ep=epoch_model)

    lat_list = list()
    lon_list = list()
    day_rad_list = list()
    measured_var_list = list()
    generated_var_list = list()
    suazade_var_list = list()
    gloria_var_list = list()
    number_seasonal_sample = 0

    sum_generated = torch.zeros(1, 1, 200)
    sum_measured = torch.zeros(1, 1, 200)
    sum_gloria = torch.zeros(1, 1, 200)
    sum_suazade = torch.zeros(200)

    for sample in ds:
        year, day_rad, lat, lon, temp, psal, doxy, measured_var = sample
        day_sample = from_day_rad_to_day(day_rad=day_rad)
        if season != "all" and not dict_season[season][0] <= day_sample <= dict_season[season][
            1] and random.random() < 0.2:
            continue
        number_seasonal_sample += 1
        generated_suazade_var = get_suazade_profile(year, day_rad, lat, lon, torch.squeeze(temp), torch.squeeze(psal),
                                                    torch.squeeze(doxy), torch.squeeze(measured_var))

        suazade_var_list.append(generated_suazade_var)
        sum_suazade += generated_suazade_var

        generated_gloria_var = get_gloria_profile(year, day_rad, lat, lon, torch.squeeze(temp), torch.squeeze(psal),
                                                  torch.squeeze(doxy), torch.squeeze(measured_var))
        gloria_var_list.append(generated_gloria_var)
        sum_gloria += generated_gloria_var

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

        generated_var_list.append(generated_var)
        sum_generated += generated_var
        measured_var_list.append(measured_var)
        sum_measured += measured_var

        depth = np.linspace(0, dict_max_pressure["NITRATE"], len(torch.squeeze(generated_var)))

        cut_sup = 3
        cut_inf = 10

        plt.plot(torch.squeeze(measured_var)[cut_sup:-cut_inf], depth[cut_sup:-cut_inf], label="measured")
        plt.plot(generated_suazade_var[cut_sup:-cut_inf], depth[cut_sup:-cut_inf], label="MLP Fourrier")
        plt.plot(torch.squeeze(generated_gloria_var)[cut_sup:-cut_inf], depth[cut_sup:-cut_inf],
                 label="MLP Pietropolli")
        plt.plot(torch.squeeze(generated_var)[cut_sup:-cut_inf], depth[cut_sup:-cut_inf], label="convolutional")
        plt.gca().invert_yaxis()

        plt.xlabel(f"{variable} {dict_unit_measure[variable]}")
        plt.ylabel('depth')

        plt.legend()
        plt.savefig(f"{path_analysis}/method_comparison_{round(lat.item(), 2)}_{round(lon.item(), 2)}.png")
        plt.close()

    generated_mean = sum_generated / number_seasonal_samples
    measured_mean = sum_measured / number_seasonal_samples
    gloria_mean = sum_gloria / number_seasonal_samples
    suazade_mean = sum_suazade / number_seasonal_samples

    max_pres = dict_max_pressure[variable]
    depth = np.linspace(0, max_pres, len(generated_profile.detach().numpy()))
    plt.plot(torch.squeeze(measured_mean), depth, label=f"measured")
    plt.plot(suazade_mean, depth, label=f"MLP Fourrier")
    plt.plot(torch.squeeze(gloria_mean), depth, label=f"MLP Pietropolli")
    plt.plot(torch.squeeze(generated_mean), depth, label=f"convolutional")
    plt.gca().invert_yaxis()

    plt.xlabel(f"{variable} [{dict_unit_measure[variable]}]")
    plt.ylabel('depth [m]')

    plt.legend()
    plt.savefig(f"{path_analysis}/method_comparison_{season}.png")
    plt.show()
    plt.close()

    return


def get_reconstruction_comparison_gpietrop(variable, season, date_model, epoch_model, mode):
    # Create savedir
    path_analysis = os.getcwd() + f"/../results/NITRATE/{date_model}/fig/comp_gpietrop"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)
    # Upload the input ds

    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}

    path_float = f"/home/gpietropolli/Desktop/canyon-float/ds/NITRATE/float_ds_sf_{mode}.csv"
    if mode == "all":
        path_float = f"/home/gpietropolli/Desktop/canyon-float/ds/NITRATE/float_ds_sf.csv"
    dataset = FloatDataset(path_float)
    ds = DataLoader(dataset, shuffle=True)

    dir_model = os.getcwd() + f"/../results/NITRATE/{date_model}/model"
    info = pd.read_csv(os.getcwd() + f"/../results/NITRATE/{date_model}/info.csv")

    # Upload and evaluate the model
    model_day, model_year, model_lat, model_lon, model = upload_and_evaluate_model(dir_model=dir_model, info_model=info,
                                                                                   ep=epoch_model)

    lat_list = list()
    lon_list = list()
    day_rad_list = list()
    measured_var_list = list()
    gloria_var_list = list()
    number_seasonal_sample = 0

    sum_measured = torch.zeros(1, 1, 200)
    sum_gloria = torch.zeros(1, 1, 200)

    for sample in ds:
        year, day_rad, lat, lon, temp, psal, doxy, measured_var = sample
        day_sample = from_day_rad_to_day(day_rad=day_rad)
        if season != "all" and not dict_season[season][0] <= day_sample <= dict_season[season][
            1] and random.random() < 0.2:
            continue
        number_seasonal_sample += 1

        generated_gloria_var = get_gloria_profile(year, day_rad, lat, lon, torch.squeeze(temp), torch.squeeze(psal),
                                                  torch.squeeze(doxy), torch.squeeze(measured_var))
        gloria_var_list.append(generated_gloria_var)
        sum_gloria += generated_gloria_var

        lat_list.append(lat.item())
        lon_list.append(lon.item())
        day_rad_list.append(day_rad.item())

        measured_var_list.append(measured_var)
        sum_measured += measured_var

        depth = np.linspace(0, dict_max_pressure["NITRATE"], len(torch.squeeze(measured_var)))

        cut_sup = 3
        cut_inf = 10

        plt.plot(torch.squeeze(measured_var)[cut_sup:-cut_inf], depth[cut_sup:-cut_inf], lw=2.5, label="measured",
                 color='green')
        plt.plot(torch.squeeze(generated_gloria_var)[cut_sup:-cut_inf], depth[cut_sup:-cut_inf], lw=2.5,
                 linestyle="dashed", label="MLP Pietropolli", color='darkorange')
        plt.gca().invert_yaxis()

        plt.xlabel(f"{variable} [{dict_unit_measure[variable]}]")
        plt.ylabel('depth [m]')

        plt.legend()
        plt.savefig(f"{path_analysis}/method_comparison_{round(lat.item(), 2)}_{round(lon.item(), 2)}.png")
        plt.close()

    measured_mean = sum_measured / number_seasonal_samples
    gloria_mean = sum_gloria / number_seasonal_samples

    max_pres = dict_max_pressure[variable]
    depth = np.linspace(0, max_pres, len(generated_profile.detach().numpy()))
    plt.plot(torch.squeeze(measured_mean), depth, lw=3, label=f"measured")
    plt.plot(torch.squeeze(gloria_mean), depth, lw=3, label=f"MLP Pietropolli")
    plt.gca().invert_yaxis()

    plt.xlabel(f"{variable} [{dict_unit_measure[variable]}]")
    plt.ylabel('depth [m]')

    plt.legend()
    plt.savefig(f"{path_analysis}/method_comparison_{season}.png")
    plt.show()
    plt.close()

    return



def seasonal_bp(variable, date_model, epoch_model, mode):
    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}
    list_loss_for_season = [[] for i in range(4)]

    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    number_samples = len(generated_var_list)
    for index_season in range(4):
        season = list(dict_season.keys())[index_season]
        for index_sample in range(number_samples):
            day_sample = from_day_rad_to_day(day_rad=day_rad_list[index_sample])
            if dict_season[season][0] <= day_sample <= dict_season[season][1]:
                loss_sample = mse_loss(generated_var_list[index_sample], measured_var_list[index_sample])
                list_loss_for_season[index_season].append(loss_sample)

    sns.boxplot(data=list_loss_for_season,
                showfliers=False)
    plt.xticks(range(4), list(dict_season.keys()))
    plt.show()
    plt.close()

    return


def geographic_bp(variable, date_model, epoch_model, mode):
    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
               'SWM': [[32, 40], [-2, 9.5]],
               'TIR': [[37, 45], [9.5, 16]],
               'ION': [[30, 37], [9.5, 22]],
               'LEV': [[30, 37], [22, 36]]}
    list_loss_ga = [[] for _ in range(len(list(dict_ga.keys())))]

    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    number_samples = len(generated_var_list)
    for index_ga in range(len(list(dict_ga.keys()))):
        ga = list(dict_ga.keys())[index_ga]
        for i in range(number_samples):
            if dict_ga[ga][0][0] <= lat_list[i] <= dict_ga[ga][0][1] and dict_ga[ga][1][0] <= lon_list[i] <= \
                    dict_ga[ga][1][1]:
                loss_sample = mse_loss(generated_var_list[i], measured_var_list[i])
                list_loss_ga[index_ga].append(loss_sample)

    sns.boxplot(data=list_loss_ga,
                showfliers=False)
    plt.xticks(range(len(list(dict_ga.keys()))), list(dict_ga.keys()))
    plt.show()
    plt.close()

    return

def seasonal_profile(season, variable, date_model, epoch_model, mode):
    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}
    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    len_profile = int(generated_var_list[0].size(dim=0))

    generated_profile = torch.zeros(len_profile)
    measured_profile = torch.zeros(len_profile)
    number_samples = len(generated_var_list)
    number_seasonal_samples = 0
    for index_sample in range(number_samples):
        day_sample = from_day_rad_to_day(day_rad=day_rad_list[index_sample])
        if dict_season[season][0] <= day_sample <= dict_season[season][1]:
            number_seasonal_samples += 1
            generated_profile = generated_profile + generated_var_list[index_sample]
            measured_profile = measured_profile + measured_var_list[index_sample]

    generated_profile = generated_profile / number_seasonal_samples
    measured_profile = measured_profile / number_seasonal_samples

    max_pres = dict_max_pressure[variable]
    depth = np.linspace(0, max_pres, len(generated_profile.detach().numpy()))
    plt.plot(generated_profile.detach().numpy(), depth, label=f"generated {variable}")
    plt.plot(measured_profile.detach().numpy(), depth, label=f"measured {variable}")
    plt.gca().invert_yaxis()

    plt.legend()
    plt.show()
    plt.close()

    return


def ga_profile(ga, variable, date_model, epoch_model, mode):
    font = {'weight': 'bold',
            'size': 100}
    matplotlib.rc('font', **font)
    plt.figure(figsize=(6.5, 8))

    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/fig/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)
    path_analysis = path_analysis + "ga_prof/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
               'SWM': [[32, 40], [-2, 9.5]],
               'TYR': [[37, 45], [9.5, 16]],
               'ION': [[30, 37], [9.5, 22]],
               'LEV': [[30, 37], [22, 36]]}

    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    len_profile = int(generated_var_list[0].size(dim=0))

    generated_profile = torch.zeros(len_profile)
    measured_profile = torch.zeros(len_profile)
    number_samples = len(generated_var_list)
    number_ga_samples = 0
    for index_sample in range(number_samples):
        if dict_ga[ga][0][0] <= lat_list[index_sample] <= dict_ga[ga][0][1] and dict_ga[ga][1][0] <= lon_list[
            index_sample] <= dict_ga[ga][1][1]:
            number_ga_samples += 1

            generated_profile = generated_profile + generated_var_list[index_sample]
            measured_profile = measured_profile + measured_var_list[index_sample]

    generated_profile = generated_profile / number_ga_samples
    measured_profile = measured_profile / number_ga_samples

    max_pres = dict_max_pressure[variable]
    depth = np.linspace(0, max_pres, len(generated_profile.detach().numpy()))
    plt.plot(generated_profile.detach().numpy(), depth, lw=3,
             linestyle="dashed", color=dict_color[ga], label=f"CNN generated")
    plt.plot(measured_profile.detach().numpy(), depth,  lw=3,
             linestyle="solid", color=dict_color[ga], label=f"measured")
    plt.gca().invert_yaxis()

    plt.xlim([0, 8.5])

    plt.title(ga)
    plt.xlabel(f"{variable} [{dict_unit_measure[variable]}]")
    plt.ylabel('depth [m]')

    plt.legend()
    plt.savefig(f"{path_analysis}profile_mean_{ga}_{mode}_{epoch_model}.png")
    plt.close()

    return


def profile_error(variable, date_model, epoch_model, mode):
    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/fig/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
               'SWM': [[32, 40], [-2, 9.5]],
               'TYR': [[37, 45], [9.5, 16]],
               'ION': [[30, 37], [9.5, 22]],
               'LEV': [[30, 37], [22, 36]]}

    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}

    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)

    fig, axs = plt.subplots(2, 2)

    for ga in list(dict_ga.keys()):
        generated, measured = get_profile_list_season_ga("W", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                                         measured_var_list)
        max_pres = dict_max_pressure[variable]
        depth = np.linspace(0, max_pres, len(generated.detach().numpy()))

        axs[0, 0].plot(np.abs(measured.detach().numpy()) -
                       np.abs(generated.detach().numpy()), depth,
                       label=ga, linestyle="solid", color=dict_color[ga])
        axs[0, 0].invert_yaxis()

        generated, measured = get_profile_list_season_ga("SP", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                                         measured_var_list)
        axs[0, 1].plot(np.abs(measured.detach().numpy()) -
                       np.abs(generated.detach().numpy()), depth,
                       label=ga, linestyle="solid", color=dict_color[ga])
        axs[0, 1].invert_yaxis()

        generated, measured = get_profile_list_season_ga("SU", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                                         measured_var_list)
        axs[1, 0].plot(np.abs(measured.detach().numpy()) -
                       np.abs(generated.detach().numpy()), depth,
                       label=ga, linestyle="solid", color=dict_color[ga])
        axs[1, 0].invert_yaxis()

        generated, measured = get_profile_list_season_ga("A", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                                         measured_var_list)
        axs[1, 1].plot(np.abs(measured.detach().numpy()) -
                       np.abs(generated.detach().numpy()), depth,
                       label=ga, linestyle="solid", color=dict_color[ga])
        axs[1, 1].invert_yaxis()

    axs[0, 0].set_title(list(dict_season.keys())[0])
    axs[0, 1].set_title(list(dict_season.keys())[1])
    axs[1, 0].set_title(list(dict_season.keys())[2])
    axs[1, 1].set_title(list(dict_season.keys())[3])

    for ax in axs.flat:
        ax.set(xlabel=f"{variable} ({dict_unit_measure[variable]})", ylabel='depth')
        # ax.set_xticks(range(1, len(list(dict_ga.keys())) + 1), list(dict_ga.keys()))

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    # for ax in axs.flat:
    #    ax.label_outer()

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize="5")

    plt.savefig(f"{path_analysis}profile_error_{mode}_{epoch_model}.png")

    # plt.show()
    plt.close()

    return


def seasonal_ga_list(season, ga, lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list):
    dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
               'SWM': [[32, 40], [-2, 9.5]],
               'TYR': [[37, 45], [9.5, 16]],
               'ION': [[30, 37], [9.5, 22]],
               'LEV': [[30, 37], [22, 36]]}

    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}

    generated_profile = []
    measured_profile = []
    number_samples = len(generated_var_list)
    for i in range(number_samples):
        day_sample = from_day_rad_to_day(day_rad=day_rad_list[i])
        if dict_season[season][0] <= day_sample <= dict_season[season][1]:
            if dict_ga[ga][0][0] <= lat_list[i] <= dict_ga[ga][0][1] and dict_ga[ga][1][0] <= lon_list[i] <= \
                    dict_ga[ga][1][1]:
                generated_profile.append(generated_var_list[i])
                measured_profile.append(measured_var_list[i])

    return generated_profile, measured_profile


def ga_variance(ga, variable, date_model, epoch_model, mode):
    dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
               'SWM': [[32, 40], [-2, 9.5]],
               'TYR': [[37, 45], [9.5, 16]],
               'ION': [[30, 37], [9.5, 22]],
               'LEV': [[30, 37], [22, 36]]}

    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    len_profile = int(generated_var_list[0].size(dim=0))

    generated_profile = []
    measured_profile = []
    number_samples = len(generated_var_list)

    for i in range(number_samples):
        if dict_ga[ga][0][0] <= lat_list[i] <= dict_ga[ga][0][1] and dict_ga[ga][1][0] <= lon_list[i] <= dict_ga[ga][1][1]:
            generated_profile.append(generated_var_list[i])
            measured_profile.append(measured_var_list[i])

    std_reconstruction = torch.zeros(len_profile)
    std_measured = torch.zeros(len_profile)

    for k in range(len_profile):
        reconstruction_depth = []
        measured_depth = []
        for prof in generated_profile:
            reconstruction_depth.append(prof[k].item())
        std_reconstruction[k] = np.std(reconstruction_depth)

        for prof in measured_profile:
            measured_depth.append(prof[k].item())
        std_measured[k] = np.std(measured_depth)

    max_pres = dict_max_pressure[variable]
    depth = np.linspace(0, max_pres, len(std_reconstruction.detach().numpy()))
    plt.plot(std_reconstruction.detach().numpy(), depth,
             linestyle="dashed", color=dict_color[ga], label=f"CNN generated")
    plt.plot(std_measured.detach().numpy(), depth,
             linestyle="solid", color=dict_color[ga], label=f"measured")

    plt.gca().invert_yaxis()

    plt.title(ga)
    plt.xlabel(f"{variable} [{dict_unit_measure[variable]}]")
    plt.ylabel('depth [m]')

    plt.legend()
    # plt.savefig(f"{path_analysis}profile_std_{ga}_{mode}_{epoch_model}.png")
    plt.show()
    plt.close()

    return


def profile_efficency(variable, date_model, epoch_model, mode):
    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/fig/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
               'SWM': [[32, 40], [-2, 9.5]],
               'TYR': [[37, 45], [9.5, 16]],
               'ION': [[30, 37], [9.5, 22]],
               'LEV': [[30, 37], [22, 36]]}

    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}

    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    fig, axs = plt.subplots(2, 2)

    for ga in list(dict_ga.keys()):
        variance = get_variance_list_season_ga("W", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                               measured_var_list)
        generated_var_list_s_ga, measured_var_list_s_ga = seasonal_ga_list("W", ga, lat_list, lon_list, day_rad_list,
                                                                           generated_var_list, measured_var_list)
        efficency = compute_efficency(generated_var_list_s_ga, measured_var_list_s_ga, variance)

        max_pres = dict_max_pressure[variable]
        depth = np.linspace(0, max_pres, len(variance.detach().numpy()))

        axs[0, 0].plot(efficency.detach().numpy(), depth, color=dict_color[ga], label=ga)
        axs[0, 0].invert_yaxis()

        variance = get_variance_list_season_ga("SP", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                               measured_var_list)
        generated_var_list_s_ga, measured_var_list_s_ga = seasonal_ga_list("SP", ga, lat_list, lon_list, day_rad_list,
                                                                           generated_var_list, measured_var_list)
        efficency = compute_efficency(generated_var_list_s_ga, measured_var_list_s_ga, variance)

        axs[0, 1].plot(efficency.detach().numpy(), depth, color=dict_color[ga], label=ga)
        axs[0, 1].invert_yaxis()

        variance = get_variance_list_season_ga("SU", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                               measured_var_list)
        generated_var_list_s_ga, measured_var_list_s_ga = seasonal_ga_list("SU", ga, lat_list, lon_list, day_rad_list,
                                                                           generated_var_list, measured_var_list)
        efficency = compute_efficency(generated_var_list_s_ga, measured_var_list_s_ga, variance)

        axs[1, 0].plot(efficency.detach().numpy(), depth, color=dict_color[ga], label=ga)
        axs[1, 0].invert_yaxis()

        variance = get_variance_list_season_ga("A", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                               measured_var_list)
        generated_var_list_s_ga, measured_var_list_s_ga = seasonal_ga_list("A", ga, lat_list, lon_list, day_rad_list,
                                                                           generated_var_list, measured_var_list)
        efficency = compute_efficency(generated_var_list_s_ga, measured_var_list_s_ga, variance)

        axs[1, 1].plot(efficency.detach().numpy(), depth, color=dict_color[ga], label=ga)
        axs[1, 1].invert_yaxis()

    axs[0, 0].set_title(list(dict_season.keys())[0])
    axs[0, 1].set_title(list(dict_season.keys())[1])
    axs[1, 0].set_title(list(dict_season.keys())[2])
    axs[1, 1].set_title(list(dict_season.keys())[3])

    for ax in axs.flat:
        ax.set(xlabel=f"{variable} ({dict_unit_measure[variable]})", ylabel='depth')
        # ax.set_xticks(range(1, len(list(dict_ga.keys())) + 1), list(dict_ga.keys()))

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize="5")

    # plt.savefig(f"{path_analysis}variance_comparison_{mode}_{epoch_model}.png")

    plt.show()
    plt.close()

    return

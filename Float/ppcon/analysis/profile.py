from utils_analysis import *


def get_profile_list_season_ga(season, ga, lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list):
    dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
               'SWM': [[32, 40], [-2, 9.5]],
               'TYR': [[37, 45], [9.5, 16]],
               'ION': [[30, 37], [9.5, 22]],
               'LEV': [[30, 37], [22, 36]]}

    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}

    len_profile = int(generated_var_list[0].size(dim=0))

    generated_profile = torch.zeros(len_profile)
    measured_profile = torch.zeros(len_profile)
    number_samples = len(generated_var_list)
    number_seasonal_samples = 0
    for i in range(number_samples):
        day_sample = from_day_rad_to_day(day_rad=day_rad_list[i])
        if dict_season[season][0] <= day_sample <= dict_season[season][1]:
            if dict_ga[ga][0][0] <= lat_list[i] <= dict_ga[ga][0][1] and dict_ga[ga][1][0] <= lon_list[i] <= \
                    dict_ga[ga][1][1]:
                number_seasonal_samples += 1
                generated_profile = generated_profile + generated_var_list[i]
                measured_profile = measured_profile + measured_var_list[i]

    generated_profile = generated_profile / number_seasonal_samples
    measured_profile = measured_profile / number_seasonal_samples

    return generated_profile, measured_profile


def profile_season_ga(variable, date_model, epoch_model, mode):
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

    fig, axs = plt.subplots(2, 2, figsize=(7, 7))

    for ga in list(dict_ga.keys()):
        generated, measured = get_profile_list_season_ga("W", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                                         measured_var_list)
        max_pres = dict_max_pressure[variable]
        depth = np.linspace(0, max_pres, len(generated.detach().numpy()))

        axs[0, 0].plot(generated.detach().numpy(), depth,
                       linestyle="dashed", color=dict_color[ga])
        axs[0, 0].plot(measured.detach().numpy(), depth,
                       label=ga, linestyle="solid", color=dict_color[ga])
        axs[0, 0].invert_yaxis()

        generated, measured = get_profile_list_season_ga("SP", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                                         measured_var_list)
        axs[0, 1].plot(generated.detach().numpy(), depth,
                       linestyle="dashed", color=dict_color[ga])
        axs[0, 1].plot(measured.detach().numpy(), depth,
                       label=ga, linestyle="solid", color=dict_color[ga])
        axs[0, 1].invert_yaxis()

        generated, measured = get_profile_list_season_ga("SU", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                                         measured_var_list)
        axs[1, 0].plot(generated.detach().numpy(), depth,
                       linestyle="dashed", color=dict_color[ga])
        axs[1, 0].plot(measured.detach().numpy(), depth,
                       label=ga, linestyle="solid", color=dict_color[ga])
        axs[1, 0].invert_yaxis()

        generated, measured = get_profile_list_season_ga("A", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                                         measured_var_list)
        axs[1, 1].plot(generated.detach().numpy(), depth,
                       linestyle="dashed", color=dict_color[ga])
        axs[1, 1].plot(measured.detach().numpy(), depth,
                       label=ga, linestyle="solid", color=dict_color[ga])
        axs[1, 1].invert_yaxis()

    axs[0, 0].set_title("Winter")
    axs[0, 1].set_title("Spring")
    axs[1, 0].set_title("Summer")
    axs[1, 1].set_title("Autumn")

    for ax in axs.flat:
        ax.set(xlabel=f"{dict_var_name[variable]} ({dict_unit_measure[variable]})", ylabel=r"Depth [$m$]")
        if variable == "BBP700":
            ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0e}'))
        # ax.set_xticks(range(1, len(list(dict_ga.keys())) + 1), list(dict_ga.keys()))

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize="5")
    plt.tight_layout()

    plt.savefig(f"{path_analysis}profile_comparison_{mode}_{epoch_model}.png")
    plt.savefig(os.getcwd() + f"/../results/paper_fig/profile_comparison_{variable}_{mode}.png", dpi=1200)

    # plt.show()
    plt.close()

    return


def get_variance_list_season_ga(season, ga, lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list):
    dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
               'SWM': [[32, 40], [-2, 9.5]],
               'TYR': [[37, 45], [9.5, 16]],
               'ION': [[30, 37], [9.5, 22]],
               'LEV': [[30, 37], [22, 36]]}

    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}

    len_profile = int(generated_var_list[0].size(dim=0))

    generated_profile = []
    # measured_profile = []
    number_samples = len(generated_var_list)
    number_seasonal_samples = 0
    for i in range(number_samples):
        day_sample = from_day_rad_to_day(day_rad=day_rad_list[i])
        if dict_season[season][0] <= day_sample <= dict_season[season][1]:
            if dict_ga[ga][0][0] <= lat_list[i] <= dict_ga[ga][0][1] and dict_ga[ga][1][0] <= lon_list[i] <= \
                    dict_ga[ga][1][1]:
                number_seasonal_samples += 1
                generated_profile.append(generated_var_list[i])
                # measured_profile.append(measured_var_list[i])

    std_reconstruction = torch.zeros(len_profile)
    for k in range(len_profile):
        reconstruction_depth = []
        for prof in generated_profile:
            reconstruction_depth.append(prof[k].item())
        # print(reconstruction_depth)
        std_reconstruction[k] = np.std(reconstruction_depth)

    return std_reconstruction


def profile_variance(variable, date_model, epoch_model, mode):
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
        max_pres = dict_max_pressure[variable]
        depth = np.linspace(0, max_pres, len(variance.detach().numpy()))

        axs[0, 0].plot(variance.detach().numpy(), depth, color=dict_color[ga], label=ga)
        axs[0, 0].invert_yaxis()

        variance = get_variance_list_season_ga("SP", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                               measured_var_list)
        axs[0, 1].plot(variance.detach().numpy(), depth, color=dict_color[ga], label=ga)
        axs[0, 1].invert_yaxis()

        variance = get_variance_list_season_ga("SU", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                               measured_var_list)
        axs[1, 0].plot(variance.detach().numpy(), depth, color=dict_color[ga], label=ga)
        axs[1, 0].invert_yaxis()

        variance = get_variance_list_season_ga("A", ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                               measured_var_list)
        axs[1, 1].plot(variance.detach().numpy(), depth, color=dict_color[ga], label=ga)
        axs[1, 1].invert_yaxis()

    axs[0, 0].set_title(list(dict_season.keys())[0])
    axs[0, 1].set_title(list(dict_season.keys())[1])
    axs[1, 0].set_title(list(dict_season.keys())[2])
    axs[1, 1].set_title(list(dict_season.keys())[3])

    for ax in axs.flat:
        ax.set(xlabel=f"{variable} ({dict_unit_measure[variable]})", ylabel=r"depth [$m$]")
        # ax.set_xticks(range(1, len(list(dict_ga.keys())) + 1), list(dict_ga.keys()))

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize="5")

    plt.savefig(f"{path_analysis}variance_comparison_{mode}_{epoch_model}.png")

    # plt.show()
    plt.close()

    return
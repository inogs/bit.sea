from utils_analysis import *


def bp_season_ga(variable, date_model, epoch_model, mode):
    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/fig/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
               'SWM': [[32, 40], [-2, 9.5]],
               'TYR': [[37, 45], [9.5, 16]],
               'ION': [[30, 37], [9.5, 22]],
               'LEV': [[30, 37], [22, 36]]}

    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}

    list_loss = [[[] for _ in range(len(list(dict_ga.keys())))] for _ in range(len(list(dict_season.keys())))]

    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    number_samples = len(generated_var_list)
    for index_s in range(len(list(dict_season.keys()))):
        season = list(dict_season.keys())[index_s]
        for index_ga in range(len(list(dict_ga.keys()))):
            ga = list(dict_ga.keys())[index_ga]
            for i in range(number_samples):
                day_sample = from_day_rad_to_day(day_rad=day_rad_list[i])
                if dict_season[season][0] <= day_sample <= dict_season[season][1]:
                    if dict_ga[ga][0][0] <= lat_list[i] <= dict_ga[ga][0][1] and dict_ga[ga][1][0] <= lon_list[i] <= \
                            dict_ga[ga][1][1]:
                        loss_sample = mse_loss(generated_var_list[i], measured_var_list[i])
                        list_loss[index_s][index_ga].append(loss_sample)

    fig, axs = plt.subplots(2, 2)

    sns.boxplot(list_loss[0],
                palette="magma",
                ax=axs[0, 0])
    axs[0, 0].set_title(list(dict_season.keys())[0])

    sns.boxplot(list_loss[1],
                palette="magma",
                ax=axs[0, 1])
    axs[0, 1].set_title(list(dict_season.keys())[1])

    sns.boxplot(list_loss[2],
                palette="magma",
                ax=axs[1, 0])
    axs[1, 0].set_title(list(dict_season.keys())[2])

    sns.boxplot(list_loss[3],
                palette="magma",
                ax=axs[1, 1])
    axs[1, 1].set_title(list(dict_season.keys())[3])

    for ax in axs.flat:
        ax.set(xlabel=f"{variable} ({dict_unit_measure[variable]})", ylabel='fitness')
        if variable == "NITRATE":
            ax.set_ylim([0.0, 1])
        if variable == "CHLA":
            ax.set_ylim([0.0, 0.08])
        if variable == "BBP700":
            ax.set_ylim([0.0, 0.00000015])
        ax.set_xticks(range(0, len(list(dict_ga.keys()))), list(dict_ga.keys()))

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

    plt.savefig(f"{path_analysis}bp_comparison_{mode}_{epoch_model}.png")

    # plt.show()
    plt.close()

    return

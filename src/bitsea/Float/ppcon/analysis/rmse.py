import matplotlib

from utils_analysis import *


def rmse(variable, date_model, epoch_model, mode):
    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    my_loss = 0
    number_samples = len(generated_var_list)
    for index_sample in range(number_samples):
        loss_sample = np.sqrt(mse_loss(generated_var_list[index_sample], measured_var_list[index_sample]))
        my_loss += loss_sample
    my_loss = my_loss / number_samples
    return my_loss


def seasonal_rmse(season, variable, date_model, epoch_model, mode):
    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}
    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    my_loss = 0
    number_samples = len(generated_var_list)
    number_seasonal_samples = 0
    for index_sample in range(number_samples):
        day_sample = from_day_rad_to_day(day_rad=day_rad_list[index_sample])
        if dict_season[season][0] <= day_sample <= dict_season[season][1]:
            number_seasonal_samples += 1
            loss_sample = np.sqrt(mse_loss(generated_var_list[index_sample], measured_var_list[index_sample]))
            my_loss += loss_sample
    my_loss = my_loss / number_seasonal_samples
    return my_loss, number_seasonal_samples


def geographic_rmse(lat_limits, lon_limits, variable, date_model, epoch_model, mode):
    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    my_loss = 0
    number_samples = len(generated_var_list)
    number_geographic_samples = 0
    for i_sample in range(number_samples):
        if lat_limits[0] <= lat_list[i_sample] <= lat_limits[1] and lon_limits[0] <= lon_list[i_sample] <= lon_limits[
            1]:
            number_geographic_samples += 1
            loss_sample = mse_loss(generated_var_list[i_sample], measured_var_list[i_sample])
            my_loss += loss_sample
    my_loss = my_loss / number_geographic_samples
    return my_loss, number_geographic_samples


def seasonal_and_geographic_rmse(variable, date_model, epoch_model, mode, make_fig=False):
    font = {'family': 'normal',
            'weight': 'bold',
            'size': 22}

    matplotlib.rc('font', **font)

    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
               'SWM': [[32, 40], [-2, 9.5]],
               'TIR': [[37, 45], [9.5, 16]],
               'ION': [[30, 37], [9.5, 22]],
               'LEV': [[30, 37], [22, 36]]}

    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}

    list_loss = [[0 for _ in range(len(list(dict_ga.keys())))] for _ in range(len(list(dict_season.keys())))]
    list_number_samples = [[0 for _ in range(len(list(dict_ga.keys())))] for _ in range(len(list(dict_season.keys())))]

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
                        loss_sample = np.sqrt(mse_loss(generated_var_list[i], measured_var_list[i]))
                        # print(list_loss)
                        list_loss[index_s][index_ga] += loss_sample
                        list_number_samples[index_s][index_ga] += 1

    list_loss = np.array(list_loss) / np.array(list_number_samples)

    # for j in range(len(list(dict_ga.keys()))):
    #     print(f"{list(dict_ga.keys())[j]} loss")
    #     loss_s = 0
    #     number_s = 0
    #     for k in range(len(list(dict_season.keys()))):
    #         loss_s += list_loss[j, k] * list_number_samples[j][k]
    #         number_s += list_number_samples[j][k]
    #     print(loss_s / number_s)

    # for j in range(len(list(dict_season.keys()))):
    #    print(f"{list(dict_season.keys())[j]} loss")
    #    loss_s = 0
    #    number_s = 0
    #    for k in range(len(list(dict_ga.keys()))):
    #        loss_s += list_loss[j, k] * list_number_samples[j][k]
    #        number_s += list_number_samples[j][k]
    #    print(loss_s / number_s)

    for j in range(len(list(dict_ga.keys()))):
        print(f"{list(dict_ga.keys())[j]} loss")
        loss_s = 0
        number_s = 0
        for k in range(len(list(dict_season.keys()))):
            loss_s += list_loss[k, j] * list_number_samples[k][j]
            number_s += list_number_samples[k][j]
        print(loss_s / number_s)

    # for j in range(len(list(dict_season.keys()))):
    # print(f"{list(dict_season.keys())[j]} losses: \t number samples")
    # for k in range(len(list(dict_ga.keys()))):
    # print(f"{list(dict_ga.keys())[k]} : {list_loss[j, k]} \t {list_number_samples[j][k]}")

    if make_fig:
        for k in range(len(dict_season.keys())):
            season = list(dict_season.keys())[k]
            plt.bar(x=range(len(dict_color.keys())), height=list_number_samples[k], color=dict_color.values())
            plt.xticks(range(len(dict_color.keys())), dict_color.keys())
            plt.ylim([0, 275])
            plt.title(f"{season} train data distribution -- {variable}")
            plt.xlabel("geographical area")
            plt.ylabel("number of samples")

            plt.savefig(f"{path_analysis}hist_{variable}_{mode}_{season}_{epoch_model}.png")
            plt.close()

    return

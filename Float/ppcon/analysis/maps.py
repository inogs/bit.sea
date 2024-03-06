from utils_analysis import *


def plot_ga_med(lat_values, lon_values, clusters, lat_limits_l=30, lat_limits_u=47, lon_limits_l=-3, lon_limits_u=37):
    lon_limits = [lon_limits_l, lon_limits_u]
    lat_limits = [lat_limits_l, lat_limits_u]
    fig = plt.figure(figsize=(8, 8))
    mediterranean_map = Basemap(llcrnrlon=lon_limits[0],
                                llcrnrlat=lat_limits[0],
                                urcrnrlon=lon_limits[1],
                                urcrnrlat=lat_limits[1],
                                resolution='h')
    mediterranean_map.drawmapboundary(fill_color='aqua')
    mediterranean_map.fillcontinents(color='coral', lake_color='aqua')
    mediterranean_map.drawcoastlines()

    plt.scatter(lon_values, lat_values, c=clusters, s=2)
    plt.ylim(lat_limits)
    plt.xlim(lon_limits)

    plt.show()
    plt.close()


def plot_med(variable, date_model, epoch_model, mode):
    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
               'SWM': [[32, 40], [-2, 9.5]],
               'TIR': [[37, 45], [9.5, 15]],
               'ION': [[30, 45], [14, 22]],
               'LEV': [[30, 37], [22, 36]]}

    dict_c = {'NWM': 0, 'SWM': 1, 'TIR': 2, 'ION': 3, 'LEV': 4}

    lat_list_season = []
    lon_list_season = []
    c_list_season = []

    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    number_samples = len(generated_var_list)

    for i in range(number_samples):
        for j in range(len(list(dict_ga.keys()))):
            key_ga = list(dict_ga.keys())[j]
            ga = dict_ga[key_ga]
            if ga[0][0] <= lat_list[i] <= ga[0][1] and ga[1][0] <= lon_list[i] <= ga[1][1]:
                lat_list_season.append(lat_list[i])
                lon_list_season.append(lon_list[i])
                c_list_season.append(dict_c[key_ga])

    plot_ga_med(lat_list_season, lon_list_season, c_list_season)


def plot_seasonal_ga_med(season, variable, date_model, epoch_model, mode):
    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}

    dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
               'SWM': [[32, 40], [-2, 9.5]],
               'TIR': [[37, 45], [9.5, 16]],
               'ION': [[30, 37], [9.5, 22]],
               'LEV': [[30, 37], [22, 36]]}

    dict_c = {'NWM': 0, 'SWM': 1, 'TIR': 2, 'ION': 3, 'LEV': 4}

    lat_list_season = []
    lon_list_season = []
    c_list_season = []

    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    number_samples = len(generated_var_list)

    for i in range(number_samples):
        day_sample = from_day_rad_to_day(day_rad=day_rad_list[i])
        if dict_season[season][0] <= day_sample <= dict_season[season][1]:

            for j in range(len(list(dict_ga.keys()))):
                key_ga = list(dict_ga.keys())[j]
                ga = dict_ga[key_ga]
                if ga[0][0] <= lat_list[i] <= ga[0][1] and ga[1][0] <= lon_list[i] <= ga[1][1]:
                    lat_list_season.append(lat_list[i])
                    lon_list_season.append(lon_list[i])
                    c_list_season.append(dict_c[key_ga])

    plot_ga_med(lat_list_season, lon_list_season, c_list_season)

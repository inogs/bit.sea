from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from analysis.utils_analysis import *
import matplotlib

dict_ga = {'NWM': [[40, 45], [-2, 9.5]],
           'SWM': [[32, 40], [-2, 9.5]],
           'TIR': [[37, 45], [9.5, 16]],
           'ION': [[30, 37], [9.5, 22]],
           'LEV': [[30, 37], [22, 36]]}

dict_season = {'W': [0, 91], 'SP': [92, 182], 'SU': [183, 273], 'A': [274, 365]}
sns.set_theme(context='paper', style='white', font='sans-serif', font_scale=1.5,
              color_codes=True, rc=None)

def seasonal_ga_variance(season, ga, lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list):
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

    return torch.mean(std_reconstruction)


def seasonal_and_geographic_std(variable, date_model, epoch_model, mode):
    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    list_std = [[0 for _ in range(len(list(dict_ga.keys())))] for _ in range(len(list(dict_season.keys())))]

    lat_list, lon_list, day_rad_list, generated_var_list, measured_var_list = get_reconstruction(variable, date_model,
                                                                                                 epoch_model, mode)
    for index_s in range(len(list(dict_season.keys()))):
        season = list(dict_season.keys())[index_s]
        for index_ga in range(len(list(dict_ga.keys()))):
            ga = list(dict_ga.keys())[index_ga]
            std_sample = seasonal_ga_variance(season, ga, lat_list, lon_list, day_rad_list, generated_var_list,
                                              measured_var_list)
            list_std[index_s][index_ga] += std_sample

    # print(list_std)
    return list_std


def seasonal_and_geographic_rmse(variable, date_model, epoch_model, mode):
    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

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
                        list_loss[index_s][index_ga] += loss_sample
                        list_number_samples[index_s][index_ga] += 1

    list_loss = np.array(list_loss) / np.array(list_number_samples)

    return list_loss, list_number_samples


def make_dim_scatter(a_list, variable):
    if variable == "NITRATE":
        list_scatter_dim = [150 if loss < 0.4 else 400 if 0.4 < loss < 0.6 else 800 for loss in a_list]
    if variable == "CHLA":
        list_scatter_dim = [150 if loss < 0.055 else 400 if 0.055 < loss < 0.1 else 800 for loss in a_list]
    if variable == "BBP700":
        list_scatter_dim = [150 if loss < 0.00022 else 400 if 0.00022 < loss < 0.0003 else 800 for loss in
                            a_list]

    return list_scatter_dim


def export_legend(legend, filename="legend.png"):
    fig = legend.figure
    fig.canvas.draw()
    bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(os.getcwd() + f"/../results/paper_fig/{filename}.png", dpi=1200, bbox_inches=bbox)


def plot_scatter_paper(variable, date_model, epoch_model, mode):

    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/fig/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    palette = ["white", "gray", "orchid", "cornflowerblue", "rebeccapurple"]
    palette = [matplotlib.colors.to_rgb(c) for c in palette]
    dict_color_blue = {'NWM': palette[0], 'SWM': palette[1], 'TYR': palette[2], 'ION': palette[3], 'LEV': palette[4]}
    # dict_color_blue = dict_color

    fig, ax = plt.subplots()
    # fig = plt.figure(figsize=(8, 5))

    list_loss, list_number_samples = seasonal_and_geographic_rmse(variable, date_model, epoch_model, "train")
    list_std = seasonal_and_geographic_std(variable, date_model, epoch_model, "train")

    patterns = ('-', '|', '.', '+')

    ax.scatter(list_std[0], list_number_samples[0],
               s=make_dim_scatter(list_loss[0], variable),
               edgecolor='black', linewidth=1,
               facecolor=[list(dict_color_blue.values())[i] + (0.55,) for i in range(len(list(dict_color_blue.values())))],
               marker="o", label="winter", hatch=2*patterns[0])  # winter
    ax.scatter(list_std[1], list_number_samples[1],
               s=make_dim_scatter(list_loss[1], variable),
               edgecolor='black', linewidth=1,
               facecolor=[list(dict_color_blue.values())[i] + (0.55,) for i in range(len(list(dict_color_blue.values())))],
               marker="o", label="spring", hatch=2*patterns[1])  # spring
    ax.scatter(list_std[2], list_number_samples[2],
               s=make_dim_scatter(list_loss[2], variable),
               edgecolor='black', linewidth=1,
               facecolor=[list(dict_color_blue.values())[i] + (0.55,) for i in range(len(list(dict_color_blue.values())))],
               marker="o", label="summer", hatch=3*patterns[2])  # summer
    ax.scatter(list_std[3], list_number_samples[3],
               s=make_dim_scatter(list_loss[3], variable),
               edgecolor='black', linewidth=1,
               facecolor=[list(dict_color_blue.values())[i] + (0.55,) for i in range(len(list(dict_color_blue.values())))],
               marker="o", label="autumn")  # autumn

    # legend 1 -- MAE dimension
    # handles, labels = scatter1.legend_elements(prop="sizes", alpha=0.6)
    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=6),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=11),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=16)
                       ]

    if variable == "NITRATE":
        un_meas = dict_unit_measure["NITRATE"]
        lg1 = ax.legend(legend_elements, ["RMSE<" + r"$0.4$", r"$0.4<$" + "RMSE" + r"$<0.6$", "RMSE" + r"$>0.6$"],
                        fontsize="10", title=f"RMSE [{un_meas}]", loc="lower right")
    if variable == "CHLA":
        un_meas = dict_unit_measure["CHLA"]
        lg1 = ax.legend(legend_elements, ["RMSE<" + r"$0.055$", r"$0.055<$" + "RMSE" + r"$<0.1$", "RMSE" + r"$>0.1$"],
                        fontsize="10", title=f"RMSE [{un_meas}]", loc="lower right")
    if variable == "BBP700":
        un_meas = dict_unit_measure["BBP700"]
        lg1 = ax.legend(legend_elements, ["RMSE<" + r"$2.2e^{-4}$", r"$2.2e^{-4}<$" + "RMSE" + r"$<3e^{-4}$",
                                          "RMSE" + r"$>3e^{-4}$"],
                        fontsize="8", title=f"RMSE [{un_meas}]", loc="lower right")

    # legend 2 -- season
    legend_elements2 = [
        Patch(facecolor="white", hatch=3*patterns[0], edgecolor='k', label="Winter"),
        Patch(facecolor="white", hatch=3*patterns[1], edgecolor='k', label="Spring"),
        Patch(facecolor="white", hatch=3*patterns[2], edgecolor='k', label="Summer"),
        Patch(facecolor="white", edgecolor='k', label="Autumn"),
    ]
    lg2 = ax.legend(handles=legend_elements2, bbox_to_anchor=(1.0, 0.35), loc="upper left", title="Season")

    # legend 3 -- ga
    legend_elements3 = [
            Patch(facecolor=list(dict_color_blue.values())[0], edgecolor='k', label=list(dict_color_blue.keys())[0]),
            Patch(facecolor=list(dict_color_blue.values())[1], edgecolor='k', label=list(dict_color_blue.keys())[1]),
            Patch(facecolor=list(dict_color_blue.values())[2], edgecolor='k', label=list(dict_color_blue.keys())[2]),
            Patch(facecolor=list(dict_color_blue.values())[3], edgecolor='k', label=list(dict_color_blue.keys())[3]),
            Patch(facecolor=list(dict_color_blue.values())[4], edgecolor='k', label=list(dict_color_blue.keys())[4]),
            ]
    # lg3 = ax.legend(handles=legend_elements3, bbox_to_anchor=(1.0, 0.35), loc="upper left", title="Geographic Area")

    ax.add_artist(lg1)
    export_legend(lg2, "legend2")

    # ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f'))
    if variable == "BBP700" or variable == "CHLA":
        ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0e}'))
    plt.xlabel(f"Standard deviation [{un_meas}]")
    plt.ylabel("Training samples")

    plt.title(dict_var_name[variable])

    plt.tight_layout()

    plt.savefig(f"{path_analysis}scatter_{mode}_{epoch_model}.png")
    plt.savefig(os.getcwd() + f"/../results/paper_fig/scatter_{variable}_{mode}.png")  #, dpi=1200)

    # plt.show()
    plt.close()

    return


def plot_scatter_paper_log(variable, date_model, epoch_model, mode):


    path_analysis = os.getcwd() + f"/../results/{variable}/{date_model}/fig/"
    if not os.path.exists(path_analysis):
        os.mkdir(path_analysis)

    palette = ["white", "gray", "orchid", "cornflowerblue", "rebeccapurple"]
    palette = [matplotlib.colors.to_rgb(c) for c in palette]
    dict_color_blue = {'NWM': palette[0], 'SWM': palette[1], 'TYR': palette[2], 'ION': palette[3], 'LEV': palette[4]}
    # dict_color_blue = dict_color

    fig, ax = plt.subplots()
    # fig = plt.figure(figsize=(8, 5))

    list_loss, list_number_samples = seasonal_and_geographic_rmse(variable, date_model, epoch_model, "train")
    list_std = seasonal_and_geographic_std(variable, date_model, epoch_model, "train")

    patterns = ('-', '|', '.', '+')

    ax.scatter(list_std[0], list_number_samples[0],
               s=make_dim_scatter(list_loss[0], variable),
               edgecolor='black', linewidth=1,
               facecolor=[list(dict_color_blue.values())[i] + (0.55,) for i in range(len(list(dict_color_blue.values())))],
               marker="o", label="winter", hatch=2*patterns[0])  # winter
    ax.scatter(list_std[1], list_number_samples[1],
               s=make_dim_scatter(list_loss[1], variable),
               edgecolor='black', linewidth=1,
               facecolor=[list(dict_color_blue.values())[i] + (0.55,) for i in range(len(list(dict_color_blue.values())))],
               marker="o", label="spring", hatch=2*patterns[1])  # spring
    ax.scatter(list_std[2], list_number_samples[2],
               s=make_dim_scatter(list_loss[2], variable),
               edgecolor='black', linewidth=1,
               facecolor=[list(dict_color_blue.values())[i] + (0.55,) for i in range(len(list(dict_color_blue.values())))],
               marker="o", label="summer", hatch=3*patterns[2])  # summer
    ax.scatter(list_std[3], list_number_samples[3],
               s=make_dim_scatter(list_loss[3], variable),
               edgecolor='black', linewidth=1,
               facecolor=[list(dict_color_blue.values())[i] + (0.55,) for i in range(len(list(dict_color_blue.values())))],
               marker="o", label="autumn")  # autumn

    ax.set_xscale('log')

    # legend 1 -- MAE dimension
    # handles, labels = scatter1.legend_elements(prop="sizes", alpha=0.6)
    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=6),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=11),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=16)
                       ]

    if variable == "NITRATE":
        un_meas = dict_unit_measure["NITRATE"]
        lg1 = ax.legend(legend_elements, ["RMSE<" + r"$0.5$", r"$0.5<$" + "RMSE" + r"$<0.7$", "RMSE" + r"$>0.7$"],
                        fontsize="10", title=f"RMSE [{un_meas}]", loc="lower right")
    if variable == "CHLA":
        un_meas = dict_unit_measure["CHLA"]
        lg1 = ax.legend(legend_elements, ["RMSE<" + r"$0.055$", r"$0.055<$" + "RMSE" + r"$<0.1$", "RMSE" + r"$>0.1$"],
                        fontsize="10", title=f"RMSE [{un_meas}]", loc="lower right")
    if variable == "BBP700":
        un_meas = dict_unit_measure["BBP700"]
        lg1 = ax.legend(legend_elements, ["RMSE<" + r"$2.2e^{-4}$", r"$2.2e^{-4}<$" + "RMSE" + r"$<3e^{-4}$",
                                          "RMSE" + r"$>3e^{-4}$"],
                        fontsize="8", title=f"RMSE [{un_meas}]", loc="lower right")

    # legend 2 -- season
    legend_elements2 = [
        Patch(facecolor="white", hatch=3*patterns[0], edgecolor='k', label="Winter"),
        Patch(facecolor="white", hatch=3*patterns[1], edgecolor='k', label="Spring"),
        Patch(facecolor="white", hatch=3*patterns[2], edgecolor='k', label="Summer"),
        Patch(facecolor="white", edgecolor='k', label="Autumn"),
    ]
    # lg2 = ax.legend(handles=legend_elements2, bbox_to_anchor=(1.0, 0.35), loc="upper left", title="Season")

    # legend 3 -- ga
    legend_elements3 = [
            Patch(facecolor=list(dict_color_blue.values())[0], edgecolor='k', label=list(dict_color_blue.keys())[0]),
            Patch(facecolor=list(dict_color_blue.values())[1], edgecolor='k', label=list(dict_color_blue.keys())[1]),
            Patch(facecolor=list(dict_color_blue.values())[2], edgecolor='k', label=list(dict_color_blue.keys())[2]),
            Patch(facecolor=list(dict_color_blue.values())[3], edgecolor='k', label=list(dict_color_blue.keys())[3]),
            Patch(facecolor=list(dict_color_blue.values())[4], edgecolor='k', label=list(dict_color_blue.keys())[4]),
            ]
    # lg3 = ax.legend(handles=legend_elements3, bbox_to_anchor=(1.0, 0.35), loc="upper left", title="Geographic Area")

    ax.add_artist(lg1)
    # ax.add_artist(lg2)
    # export_legend(lg2, "legend2")
    # ax.add_artist(lg3)
    # export_legend(lg3, "legend3")

    # ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f'))
    if variable == "BBP700" or variable == "CHLA":
        ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0e}'))
    plt.xlabel(f"Logarithmic standard deviation")
    plt.ylabel("Training samples")

    plt.title(dict_var_name[variable])
    ax.set_xticks([])

    plt.tight_layout()

    plt.savefig(f"{path_analysis}scatter_{mode}_{epoch_model}.png")
    plt.savefig(os.getcwd() + f"/../results/paper_fig/scatter_{variable}_{mode}_log.png", dpi=1200)

    # plt.show()
    plt.close()

    return

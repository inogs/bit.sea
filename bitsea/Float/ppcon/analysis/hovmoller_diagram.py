import os

import numpy as np
import netCDF4 as nc

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error

from discretization import *
from analysis.utils_analysis import dict_unit_measure
from make_ds.make_superfloat_ds import discretize

sns.set_theme(context='paper', style='white', font='sans-serif', font_scale=1.5,
              color_codes=True, rc=None)

dir_path = os.getcwd() + f"/../ds/SUPERFLOAT_PPCon/"
dir_list = os.listdir(dir_path)

for var in ["NITRATE", "CHLA", "BBP700"]:
    for folder_name in dir_list:  #dir_list:  # : dir_list
        if folder_name[0] == "F" or folder_name[0] == ".":
            continue
        folder_path = dir_path + folder_name

        files = os.listdir(folder_path)  # get ordered list of measurements
        files.sort()
        print(files)

        matrix_measured = np.zeros((200, len(files))) # initialize the matrix for the cmap function
        matrix_generated = np.zeros((200, len(files)))

        index_discarded = [] # for each of the measurements get the original measurements and the prediction
        flag_print = 0
        x_ticks = []  # date of sampling corresponding to the sampling value
        for index in range(len(files)):

            nc_file = files[index]
            nc_path = folder_path + "/" + nc_file
            try:
                ds = nc.Dataset(nc_path, "r")
            except Exception as error:
                continue

            date = ds["REFERENCE_DATE_TIME"][:]
            date = [str(np.ma.getdata(date)[k])[2] for k in range(len(date))]
            date = date[4] + date[5] + "/" + date[0] + date[1] + date[2] + date[3]

            if int(date[3] + date[4] + date[5] + date[6]) > 2020:
                index_discarded.append(index)
                continue
            else:
                x_ticks.append(date)

            if not flag_print:
                lat = float(ds["LATITUDE"][:])
                lon = float(ds["LONGITUDE"][:])

                if var == "CHLA":
                    if lon < 12:
                        maxmax = 0.75
                    else:
                        maxmax = 0.5
                    minmin = -0.02
                if var == "BBP700":
                    maxmax = 0.004
                    minmin = 0.000

                flag_print = 1

            if "DOXY" not in ds.variables.keys():
                matrix_measured[:, index] = -999 * np.array(200)
                matrix_generated[:, index] = -999 * np.array(200)
                index_discarded.append(index)
                continue

            if var not in ds.variables.keys() or f"{var}_PPCON" not in ds.variables.keys():
                matrix_generated[:, index] = -999 * np.array(200)
                matrix_measured[:, index] = -999 * np.array(200)
                index_discarded.append(index)
                continue

            if len(ds[f"PRES_{var}"][:].data) == 200:
                matrix_measured[:, index] = -999 * np.array(200)
                matrix_generated[:, index] = -999 * np.array(200)
                index_discarded.append(index)
            else:
                var_measured = ds[var][:].data
                pres_var_measured = ds[f"PRES_{var}"][:].data

                var_measured_interpolated = discretize(pres_var_measured, var_measured, dict_max_pressure[var],
                                                       dict_interval[var])

                matrix_measured[:, index] = var_measured_interpolated

                var_generated = ds[f"{var}_PPCON"][:].data
                pres_var_generated = ds[f"PRES_{var}_PPCON"][:].data

                matrix_generated[:, index] = var_generated

        try:  # find minimum of minima & maximum of maxima
            sorted_indices = sorted(index_discarded, reverse=True)
            for ind in sorted_indices:
                x_ticks.pop(ind)
        except Exception as error:
            continue


        matrix_generated_full = np.delete(matrix_generated, index_discarded, axis=1)
        matrix_measured_full = np.delete(matrix_measured, index_discarded, axis=1)

        matrix_generated_full[matrix_generated_full < 0] = 0
        matrix_measured_full[matrix_measured_full < 0] = 0

        # mse = mean_squared_error(matrix_generated_full, matrix_measured_full)
        # print(np.sqrt(mse))

        try:  # find minimum of minima & maximum of maxima
            if var == "NITRATE":
                minmin = -0.15
                max = np.max([np.max(matrix_measured), np.max(matrix_generated), 0])
                maxmax = np.min([8, max - 0.5])
        except Exception as error:
            continue

        try:
            if np.max(matrix_measured) <= 0:
                continue
        except Exception as error:
            continue

        x_number_ticks = np.arange(0, len(x_ticks))
        fig, axs = plt.subplots(2, figsize=(8.8, 6))

        cmap = matplotlib.cm.get_cmap('viridis').copy()
        cmap.set_under('white')

        im1 = axs[0].imshow(matrix_measured_full, vmin=minmin, vmax=maxmax,
                            cmap=cmap,
                            aspect='auto')  # , interpolation="nearest")
        axs[0].set_title("Measured")
        axs[0].set_xticks([])
        axs[0].set_yticks(np.arange(0, 200)[::50], np.arange(0, dict_max_pressure[var], dict_interval[var])[::50])
        axs[0].set_ylabel(r"depth [$m$]")

        im2 = axs[1].imshow(matrix_generated_full, vmin=minmin, vmax=maxmax,
                            cmap=cmap,
                            aspect='auto')  # , interpolation="bilinear")
        axs[1].set_title("PPCon prediction")
        try:
            axs[1].set_xticks(x_number_ticks[::20], x_ticks[::20], rotation=45)
        except Exception as error:
            continue
        axs[1].tick_params(axis='x', labelsize=10)
        axs[1].set_yticks(np.arange(0, 200)[::50], np.arange(0, dict_max_pressure[var], dict_interval[var])[::50])
        axs[1].set_ylabel(r"depth [$m$]")

        fig.suptitle(f"{folder_name}")
        plt.tight_layout()

        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.88, 0.15, 0.025, 0.7])
        cb = fig.colorbar(im1,  cax=cbar_ax, label=f"{var} [{dict_unit_measure[var]}]", shrink=0.1) #

        if var == "BBP700":
            cb_ticks = [float(cb.ax.get_yticklabels()[el].get_text()) for el in range(len(cb.ax.get_yticklabels()))]
            cb.ax.set_yticklabels(['{:,.0e}'.format(x) for x in cb_ticks], fontsize=8)

        if not os.path.exists(os.getcwd() + f"/../results/hovm"):
            os.mkdir(os.getcwd() + f"/../results/hovm")
        if not os.path.exists(os.getcwd() + f"/../results/hovm/{var}/"):
            os.mkdir(os.getcwd() + f"/../results/hovm/{var}/")
        plt.savefig(os.getcwd() + f"/../results/hovm/{var}/{folder_name}_{round(lat, 2)}_{round(lon, 2)}.png")  #, dpi=1200)

        plt.close()

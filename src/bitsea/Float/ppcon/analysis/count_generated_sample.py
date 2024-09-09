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

# I need as input the folder which contains float measurements and then take all the float vector one after each other

sns.set_theme(context='paper', style='white', font='sans-serif', font_scale=1.5,
              color_codes=True, rc=None)

# dir_list = os.listdir(f"/Users/admin/Desktop/ppcon/ds/SUPERFLOAT_PPCon/")
dir_path = os.getcwd() + f"/../ds/SUPERFLOAT_PPCon/"
dir_list = os.listdir(dir_path)

qc_generated = np.ones(200) * -9999

for var in ["CHLA"]:
    number_var_original = 0
    number_var_generated = 0

    for folder_name in dir_list:  # : dir_list
        if folder_name[0] == "F" or folder_name[0] == ".":
            continue
        folder_path = dir_path + folder_name

        # get ordere list of measurements
        files = os.listdir(folder_path)
        files.sort()

        for index in range(len(files)):

            nc_file = files[index]
            nc_path = folder_path + "/" + nc_file
            try:
                ds = nc.Dataset(nc_path, "r")
            except Exception as error:
                continue

            if var in ds.variables.keys():
                if f"QC_{var}" in ds.variables.keys():
                    number_var_generated += 1
                if f"{var}_QC" in ds.variables.keys():
                    number_var_original += 1

    print(var)
    print(f"number original samples {number_var_original}")
    print(f"number generated samples {number_var_generated}")
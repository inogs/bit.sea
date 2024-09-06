import os

import netCDF4 as nc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from make_ds.make_superfloat_ds import discretize
from discretization import dict_max_pressure, dict_interval

sns.set_theme(context='notebook', style='whitegrid', palette='deep', font='sans-serif', font_scale=1,
              color_codes=True, rc=None)

dict_unit_measure = {"TEMP": "°C",
                     "PSAL": "PSU",
                     "DOXY": "mg/l",
                     "NITRATE": "mmol/m^3",
                     "CHLA": "mg/m^3",
                     "BBP700": " "}

path = os.getcwd() + '/../ds/SUPERFLOAT/Float_Index.txt'
name_list = pd.read_csv(path, header=None).to_numpy()[:, 0].tolist()

i = 2900
flag = 0
var_name = "BBP700"
while flag == 0:
    path = os.getcwd() + "/../ds/SUPERFLOAT/" + name_list[i]
    # print(path)
    if not os.path.exists(path):
        i += 1
        continue

    ds = nc.Dataset(path)
    if var_name not in ds.variables.keys():
        i += 1
        continue
    else:
        flag = 1
        variables = ds.variables
        # print(variables)
        qf = ds[f"{var_name}"][:].data[:]

        # if 2 not in qf:
        #    flag = 0

        var_df = ds[var_name][:].data[:]
        pres_var = ds[f"PRES_{var_name}"][:].data[:]

# print(qf)
qf = qf
print(len(qf))
# variable = discretize(pres_var, var_df, max_pres=dict_max_pressure[var_name], interval=dict_interval[var_name])
# variable = variable * 100
# print("=====")
# print(variable)
plt.plot(qf, pres_var, linewidth=1)
plt.plot(qf, pres_var, "r.")
plt.gca().invert_yaxis()
plt.xlabel(f"{var_name} ({dict_unit_measure[var_name]})")
plt.ylabel("depth (m)")
# plt.title(var_name)
plt.savefig(os.getcwd() + f"/../results/other_fig/profile{var_name}.png")
plt.show()
plt.close()

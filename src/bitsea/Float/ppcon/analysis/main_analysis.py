from analysis.profile import profile_season_ga
# from analysis.bp import seasonal_and_geographic_bp
from analysis.rmse import *
from analysis.scatter_error import plot_scatter_paper, plot_scatter_paper_log
from analysis.comparison_architecture import *
# from maps import *
# from all_toghether import *

# reconstruction_profile_MLP("NITRATE", "all", "2023-04-04_", 50, "test")

dict_models = {
    "NITRATE": ["2023-12-16_", 100],
    "CHLA": ["2023-12-17", 150],
    "BBP700": ["2023-12-15", 125]
}

for var in ["NITRATE"]:
    date = dict_models[var][0]
    epoch = dict_models[var][1]
    reconstruction_profile_MLP(var, date, epoch, "test")

dict_comparison = {
    1: ["2023-12-15", 125],
    2: ["2023-12-17", 150],
    3: ["2023-12-17", 200]}
# my_var = "NITRATE"
# plot_med(my_var, dict_models[my_var][0], dict_models[my_var][1], "train")


for index in range(1, 4):
    var = "BBP700"
    date = dict_comparison[index][0]
    epoch = dict_comparison[index][1]
    # reconstruction_profiles(var, date, epoch, "test")
    # print(my_rmse)
    pass

for var in ["NITRATE"]:
    date = dict_models[var][0]
    epoch = dict_models[var][1]


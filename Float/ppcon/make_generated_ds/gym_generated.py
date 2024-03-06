import os

import pandas as pd
import netCDF4 as nc
from torch.utils.data import DataLoader

from dataset_with_float_names import FloatDataset
from discretization import *
import numpy as np

main_dir = "6902901"
name_file = "MR6902901_053"
path = f"/home/gpietropolli/Desktop/canyon-float/ds/SUPERFLOAT_PPCon/{main_dir}/{name_file}.nc"

ds = nc.Dataset(path)
print(list(ds.dimensions.keys()))

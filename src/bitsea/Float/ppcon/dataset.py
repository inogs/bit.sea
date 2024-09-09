import pandas as pd
import torch
from torch.utils.data import Dataset


def from_string_to_tensor(string):
    threshold = 99999
    label = 1  # label equal to one means that the sample is good
    string = string[8:-2].split(",")
    out = torch.zeros(200)
    for ind in range(len(string)):
        if float(string[ind]) >= threshold:
            label = 0
        out[ind] = torch.tensor(float(string[ind]))
    return out, label


class FloatDataset(Dataset):

    def __init__(self, path_df=None):
        super().__init__()
        if path_df is not None:
            self.path_df = path_df
            self.df = pd.read_csv(self.path_df)
        else:
            raise Exception("Paths should be given as input to initialize the Float class.")

    def __len__(self):
        """Denotes the total number of samples"""
        return len(self.df.iloc[0, :])

    def __getitem__(self, index):
        """Generates one sample of data"""
        try:
            self.samples = self.df.iloc[:, index + 1].tolist()  # Select sample
        except Exception as error:
            pass

        # self.samples = self.df.iloc[:, index + 1].tolist()  # Select sample

        year = torch.tensor(float(self.samples[0]))
        day_rad = torch.tensor(float(self.samples[1]))
        lat = torch.tensor(float(self.samples[2]))
        lon = torch.tensor(float(self.samples[3]))
        temp, label_temp = from_string_to_tensor(self.samples[4])
        psal, label_psal = from_string_to_tensor(self.samples[5])
        doxy, label_doxy = from_string_to_tensor(self.samples[6])
        nitrate, label_nitrate = from_string_to_tensor(self.samples[7])

        label = label_doxy * label_psal * label_temp  # the label is equal to one only if I have data valid for all
        # depth

        return year, day_rad, lat, lon, temp, psal, doxy, nitrate

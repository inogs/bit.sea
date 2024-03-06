import numpy as np


def preparation_dataset(data):
    for i in range(len(data[:, 0])):  # iteration on the rows (i.e. the samples)
        data[i, 1] = data[i, 1] / 90
        data[i, 2] = np.abs(1 - np.mod(data[i, 2] - 110, 360) / 180)  # fix longitude input
        data[i, 3] = np.abs(1 - np.mod(data[i, 3] - 20, 360) / 180)  # fix longitude input
        data[i, 7] = data[i, 7] / 20000 + (1 / ((1 + np.exp(-data[i, 7] / 300)) ** 3))  # fix depth input
    return data


def normalization_training(data, mean, std):
    for i in range(data.size()[1]):  # iterations over the columns
        data[:, i] = 2 / 3 * (data[:, i] - mean[i]) / std[i]
    return data
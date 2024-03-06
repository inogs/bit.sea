import torch.nn as nn


class MLPDay(nn.Module):
    def __init__(self):
        super(MLPDay, self).__init__()
        self.network = nn.Sequential(
            nn.Linear(1, 80),
            nn.SELU(),
            nn.Linear(80, 140),
            nn.SELU(),
            nn.Linear(140, 200),
            nn.SELU()
        )

    def forward(self, x):
        # x = x.view(-1, self.input_size)
        return self.network(x)


class MLPYear(nn.Module):
    def __init__(self):
        super(MLPYear, self).__init__()
        self.network = nn.Sequential(
            nn.Linear(1, 80),
            nn.SELU(),
            nn.Linear(80, 140),
            nn.SELU(),
            nn.Linear(140, 200),
            nn.SELU()
        )

    def forward(self, x):
        # x = x.view(-1, self.input_size)
        return self.network(x)


class MLPLat(nn.Module):
    def __init__(self):
        super(MLPLat, self).__init__()
        self.network = nn.Sequential(
            nn.Linear(1, 80),
            nn.SELU(),
            nn.Linear(80, 140),
            nn.SELU(),
            nn.Linear(140, 200),
            nn.SELU()
        )

    def forward(self, x):
        # x = x.view(-1, self.input_size)
        return self.network(x)


class MLPLon(nn.Module):
    def __init__(self):
        super(MLPLon, self).__init__()
        self.network = nn.Sequential(
            nn.Linear(1, 80),
            nn.SELU(),
            nn.Linear(80, 140),
            nn.SELU(),
            nn.Linear(140, 200),
            nn.SELU()
        )

    def forward(self, x):
        # x = x.view(-1, self.input_size)
        return self.network(x)
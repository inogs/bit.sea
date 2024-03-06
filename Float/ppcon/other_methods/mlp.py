"""
Implementation of the Bayesian Multi-Layer Percepton
"""
import torch
import torch.nn as nn

A = 5 / 3  # 4 / 3


def mysigmoid(x):
    return A * torch.tanh(x)  # A*(np.exp(a*x) -1)/(np.exp(a*x)+1)


class MySigmoid(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x):
        return mysigmoid(x)


F = 3


def finalsigmoid(x):
    return F * torch.tanh(x)  # A*(np.exp(a*x) -1)/(np.exp(a*x)+1)


class FinalSigmoid(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x):
        return mysigmoid(x)


activation_function = MySigmoid()
final_sigmoid = FinalSigmoid()

topology_ = [[8, 31, 23, 1], [8, 50, 30, 1], [8, 40, 20, 1], [8, 15, 8, 1], [8, 27, 23, 1], [8, 25, 25, 1],
             [8, 19, 11, 1],
             [8, 41, 9, 1], [8, 33, 29, 1], [8, 22, 15, 1]]

topology_2 = [[8, 31, 23, 1], [8, 50, 30, 1], [8, 40, 20, 1], [8, 65, 41, 1], [8, 35, 23, 1], [8, 45, 25, 1],
              [8, 29, 11, 1],
              [8, 41, 14, 1], [8, 33, 29, 1], [8, 29, 15, 1]]


class FinalSigmoid(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x):
        return mysigmoid(x)


class MLP(nn.Module):
    def __init__(self, k):
        self.topology = topology_2[k]
        input_size, n_hidden1, n_hidden2, output_size = self.topology
        super(MLP, self).__init__()
        self.input_size = input_size
        self.network = nn.Sequential(
            nn.Linear(input_size, n_hidden1),  # BayesianLinear(input_size, n_hidden1),  #
            activation_function,  # nn.SELU(),  #  funzione bene con nn.ReLU(), ELU(),
            nn.Linear(n_hidden1, n_hidden2),  # BayesianLinear(n_hidden1, n_hidden2),  #
            activation_function,  # nn.SELU(), #
            nn.Linear(n_hidden2, output_size),  # BayesianLinear(n_hidden2, output_size),  #
            nn.SELU()  # NN.SELU()
        )

    def forward(self, x):
        x = x.view(-1, self.input_size)
        return self.network(x)


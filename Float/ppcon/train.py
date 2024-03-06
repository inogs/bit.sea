import os

import numpy as np
from IPython import display

from pytorchtools import EarlyStopping

import torch
from torch.optim import Adadelta
from torch.nn.functional import mse_loss

# from conv1med import Conv1dMed
from architecture.conv1med_dp import Conv1dMed
from architecture.mlp import MLPDay, MLPYear, MLPLat, MLPLon


def train_model(train_loader, val_loader, epoch, lr, dp_rate, lambda_l2_reg, alpha_smooth_reg, attention_max, snaperiod,
                device, dir, verbose=False):

    save_dir = dir + "/model/"
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    model_mlp_day = MLPDay()
    model_mlp_year = MLPYear()
    model_mlp_lat = MLPLat()
    model_mlp_lon = MLPLon()
    model_conv = Conv1dMed(dp_rate=dp_rate)

    # Moving models to GPU when available
    model_mlp_day.to(device)
    model_mlp_year.to(device)
    model_mlp_lat.to(device)
    model_mlp_lon.to(device)
    model_conv.to(device)

    params = list(model_mlp_day.parameters()) + list(model_mlp_year.parameters()) + list(
        model_mlp_lat.parameters()) + list(
        model_mlp_lon.parameters()) + list(model_conv.parameters())
    optimizer = Adadelta(params=params, lr=lr)

    f, f_test = open(dir + "/train_loss.txt", "w+"), open(dir + "/test_loss.txt", "w+")
    f_mse_train, f_mse_test = open(dir + "/mse_train_loss.txt", "w+"), open(dir + "/mse_test_loss.txt", "w+")

    # initialize the early_stopping object
    path_checkpoint = save_dir + 'checkpoint.pt'
    early_stopping = EarlyStopping(patience=5, verbose=True, path=path_checkpoint, )

    for ep in range(epoch + 1):
        loss_train = []
        loss_test = []
        mse_train = []
        mse_test = []
        # Models in training mode
        model_mlp_day.train()
        model_mlp_year.train()
        model_mlp_lat.train()
        model_mlp_lon.train()
        model_conv.train()

        for training_year, training_day, training_lat, training_lon, training_temp, training_psal, training_doxy, training_output in train_loader:

            # Moving tensors to GPU when available
            training_year = training_year.to(device)
            training_day = training_day.to(device)
            training_lat = training_lat.to(device)
            training_lon = training_lon.to(device)
            training_temp = training_temp.to(device)
            training_psal = training_psal.to(device)
            training_doxy = training_doxy.to(device)
            training_output = training_output.to(device)

            training_day = training_day.unsqueeze(1)
            output_day = model_mlp_day(training_day)

            training_year = training_year.unsqueeze(1)
            output_year = model_mlp_year(training_year.float())

            training_lat = training_lat.unsqueeze(1)
            output_lat = model_mlp_lat(training_lat.float())

            training_lon = training_lon.unsqueeze(1)
            output_lon = model_mlp_lon(training_lon.float())

            output_day = torch.transpose(output_day.unsqueeze(0), 0, 1)
            output_year = torch.transpose(output_year.unsqueeze(0), 0, 1)
            output_lat = torch.transpose(output_lat.unsqueeze(0), 0, 1)
            output_lon = torch.transpose(output_lon.unsqueeze(0), 0, 1)
            training_temp = torch.transpose(training_temp.unsqueeze(0), 0, 1)
            training_psal = torch.transpose(training_psal.unsqueeze(0), 0, 1)
            training_doxy = torch.transpose(training_doxy.unsqueeze(0), 0, 1)
            training_output = torch.transpose(training_output.unsqueeze(0), 0, 1)

            training_x = torch.cat(
                (output_day, output_year, output_lat, output_lon, training_temp, training_psal, training_doxy), 1)

            output = model_conv(training_x.float())

            mse = mse_loss(training_output, output)

            # adding loss attention to the maximum of the peak
            max_training_output = torch.max(training_output)
            max_output = torch.max(output)
            peak_difference = attention_max * torch.abs(max_training_output - max_output)

            l2_norm = sum(p.pow(2.0).sum() for p in model_conv.parameters())
            l2_reg = lambda_l2_reg * l2_norm

            smoothness = 0
            for index in range(output.shape[0]):
                batch_tens = output[index, 0, :]
                batch_tens_smoothness = sum(torch.abs(batch_tens[i] - batch_tens[i - 1]) for i in range(1, batch_tens.shape[0]))
                smoothness += batch_tens_smoothness
            smoothness = alpha_smooth_reg * smoothness

            loss_conv = mse + l2_reg + smoothness + peak_difference

            loss_train.append(loss_conv.item())
            mse_train.append(mse.item())

            if verbose:
                print(f"[EPOCH]: {ep + 1}, [LOSS]: {loss_conv.item():.12f}")
                display.clear_output(wait=True)

            optimizer.zero_grad()
            loss_conv.backward()
            optimizer.step()

        avg_train_loss = np.average(loss_train)
        avg_train_mse = np.average(mse_train)

        print(f"[==== EPOCH]: {ep + 1}, [AVERAGE LOSS]: {avg_train_loss:.5f}")
        f.write(f"[EPOCH]: {ep + 1}, [LOSS]: {avg_train_loss:.5f} \n")
        f_mse_train.write(f"[EPOCH]: {ep + 1}, [LOSS]: {avg_train_mse:.4f} \n")

        # early_stopping needs the training loss to check if it has decreased,
        # and if it has, it will make a checkpoint of the current model
        # early_stopping(avg_train_loss, model_conv)
        # stop if validation loss doesn't improve after a given patience
        # if early_stopping.early_stop:
        #    print("Early stopping")
        #    break

        # Saving model and testing
        if ep % snaperiod == 0 or ep == epoch:

            # Saving models at epoch ep
            torch.save(model_mlp_day.state_dict(), save_dir + "/model_day_" + str(ep) + ".pt")
            torch.save(model_mlp_year.state_dict(), save_dir + "/model_year_" + str(ep) + ".pt")
            torch.save(model_mlp_lat.state_dict(), save_dir + "/model_lat_" + str(ep) + ".pt")
            torch.save(model_mlp_lon.state_dict(), save_dir + "/model_lon_" + str(ep) + ".pt")
            torch.save(model_conv.state_dict(), save_dir + "/model_conv_" + str(ep) + ".pt")

            # Model in validation mode
            model_mlp_day.eval()
            model_mlp_year.eval()
            model_mlp_lat.eval()
            model_mlp_lon.eval()
            model_conv.eval()

            with torch.no_grad():
                for testing_year, testing_day, testing_lat, testing_lon, testing_temp, testing_psal, testing_doxy, testing_output in val_loader:

                    # Moving tensors to GPU when available
                    testing_year = testing_year.to(device)
                    testing_day = testing_day.to(device)
                    testing_lat = testing_lat.to(device)
                    testing_lon = testing_lon.to(device)
                    testing_temp = testing_temp.to(device)
                    testing_psal = testing_psal.to(device)
                    testing_doxy = testing_doxy.to(device)
                    testing_output = testing_output.to(device)

                    testing_day = testing_day.unsqueeze(1)
                    output_day_test = model_mlp_day(testing_day)

                    testing_year = testing_year.unsqueeze(1)
                    output_year_test = model_mlp_year(testing_year)

                    testing_lat = testing_lat.unsqueeze(1)
                    output_lat_test = model_mlp_lat(testing_lat)

                    testing_lon = testing_lon.unsqueeze(1)
                    output_lon_test = model_mlp_lon(testing_lon)

                    output_day_test = torch.transpose(output_day_test.unsqueeze(0), 0, 1)
                    output_year_test = torch.transpose(output_year_test.unsqueeze(0), 0, 1)
                    output_lat_test = torch.transpose(output_lat_test.unsqueeze(0), 0, 1)
                    output_lon_test = torch.transpose(output_lon_test.unsqueeze(0), 0, 1)
                    testing_temp = torch.transpose(testing_temp.unsqueeze(0), 0, 1)
                    testing_psal = torch.transpose(testing_psal.unsqueeze(0), 0, 1)
                    testing_doxy = torch.transpose(testing_doxy.unsqueeze(0), 0, 1)
                    testing_output = torch.transpose(testing_output.unsqueeze(0), 0, 1)

                    testing_x = torch.cat((output_day_test, output_year_test, output_lat_test, output_lon_test,
                                           testing_temp, testing_psal, testing_doxy), 1)

                    output_test = model_conv(testing_x.float())

                    mse = mse_loss(testing_output, output_test)

                    l2_norm = sum(p.pow(2.0).sum() for p in model_conv.parameters())
                    l2_reg = lambda_l2_reg * l2_norm

                    smoothness = 0
                    for index in range(output_test.shape[0]):
                        batch_tens = output_test[index, 0, :]
                        batch_tens_smoothness = sum(batch_tens[i] - batch_tens[i-1] for i in range(1, batch_tens.shape[0]))
                        smoothness += batch_tens_smoothness
                    smoothness = alpha_smooth_reg * smoothness

                    loss_conv = mse + l2_reg + smoothness

                    loss_test.append(loss_conv)
                    mse_test.append(mse)

                    if verbose:
                        print(f"-----[EPOCH]: {ep + 1}, [TEST LOSS]: {loss_conv.item():.12f}")
                        display.clear_output(wait=True)

            avg_test_loss = np.average([loss_.cpu() for loss_ in loss_test])
            avg_test_mse = np.average([loss_.cpu() for loss_ in mse_test])

            print(f"[==== EPOCH]: {ep + 1}, [AVERAGE TEST LOSS]: {avg_test_loss:.5f}")
            f_test.write(f"[EPOCH]: {ep + 1}, [TEST LOSS]: {avg_test_loss:.5f} \n")
            f_mse_test.write(f"[EPOCH]: {ep + 1}, [TEST LOSS]: {avg_test_mse:.4f} \n")

    f.close()
    f_test.close()
    f_mse_train.close()
    f_mse_train.close()

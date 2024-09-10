import numpy as np
import torch
import datetime as dt
from datetime import datetime, date
import os, sys

from discretization import dict_max_pressure, dict_interval

basedir = os.getcwd()


def get_nitrate(timeobj, lat, lon, pres, temp, psal, doxy):
    """
    Calculates the nitrate value of canyon_med
    Arguments:
    * timeobj * a datetime object
    * lat     * scalar value
    * lon     * idem
    * pres    * idem
    * psal    * idem
    * doxy    * idem
    """
    inwgts = np.loadtxt(basedir + '/../other_methods/wgts_NO3.txt')

    lon = float(lon)
    lat = float(lat)
    pres = float(pres)
    temp = float(temp)
    psal = float(psal)
    doxy = float(doxy)
    doy = int(timeobj.strftime('%j')) * 360. / 365

    if lon > 180: lon = lon - 360
    data = [lat / 90, np.abs(1 - np.mod(lon - 110, 360) / 180.), np.abs(1 - np.mod(lon - 20, 360) / 180.), temp, psal,
            doxy, pres / 2e4 + 1. / (1 + np.exp(-pres / 300.)) ** 3]
    data = np.array(data)
    no = 1
    nol = 1

    noparsets = inwgts.shape[1] - 1
    ni = len(data)
    mw = inwgts[0:ni + 1, -1]
    sw = inwgts[ni + 1:2 * ni + 2, -1]
    data_N = (data - mw[:ni]) / sw[:ni]
    wgts = inwgts[3, 0:noparsets]
    betaciw = inwgts[2 * ni + 2:, noparsets]
    ii = np.isnan(betaciw)
    betaciw = betaciw[~ii]

    cval = np.ones((nol, noparsets), np.float32) * np.nan
    cvalcy = np.ones((1, noparsets), np.float32) * np.nan
    inval = np.ones((nol, ni, noparsets), np.float32) * np.nan

    for l in range(noparsets):
        nlayerflag = 1 + bool(inwgts[1, l])
        nl1 = int(inwgts[0, l])
        nl2 = int(inwgts[1, l])
        beta = inwgts[2, l]
        pos_1 = nl1 * ni + 4
        pos_2 = pos_1 + nl1
        pos_3 = pos_2 + nl2 * nl1
        pos_4 = pos_3 + nl2
        w1 = inwgts[4:pos_1, l].reshape(ni, nl1).T
        b1 = inwgts[pos_1:pos_2, l]
        w2 = inwgts[pos_2:pos_3, l].reshape(nl1, nl2).T
        b2 = inwgts[pos_3:pos_4, l]

        if nlayerflag == 2:
            pos_5 = pos_4 + no * nl2
            pos_6 = pos_5 + no
            w3 = inwgts[pos_4:pos_5, l].reshape(nl2, no).T
            b3 = inwgts[pos_5:pos_6, l]
        a = np.dot(data_N.T, w1.T) + b1
        if nlayerflag == 1:
            y = np.dot(np.tanh(a.T), w2.T) + b2
        if nlayerflag == 2:
            b = np.dot(np.tanh(a.T), w2.T) + b2
            y = np.dot(np.tanh(b.T), w3.T) + b3

        cval[:, l] = y
        cvalcy[:, l] = 1. / beta
        x1 = w1 * (1 - np.tanh(a) ** 2).reshape(nl1, 1)
        if nlayerflag == 1:
            inx = np.dot(w2, x1)
        if nlayerflag == 2:
            x2 = w2 * (1 - np.tanh(b) ** 2).reshape(nl2, 1)
            inx = np.dot(w3, x2).dot(x1)
        inval[:, :, l] = inx

    # Denormalization of the network output
    cval = cval * sw[ni] + mw[ni]  # variable
    cvalcy = cvalcy * sw[ni] ** 2  # 'noise' variance
    # add committee of all params as evidence-weighted mean
    V1 = wgts.sum()
    V2 = (wgts ** 2).sum()
    out_value = (wgts * cval).sum() / V1  # weighted mean
    return out_value


def get_suazade_profile(year, day_rad, lat, lon, temp, psal, doxy, measured_var):

    output_suazade = torch.zeros(len(measured_var))
    pres = np.arange(0, dict_max_pressure["NITRATE"], dict_interval["NITRATE"])
    day_total = day_rad * 365 / (2 * np.pi)
    # print(day_total % 30)
    datetime_object = dt.date(int(year.item()), int((day_total-day_total % 31)/31) + 1, 15)
    for index in range(len(measured_var)):
        output_suazade[index] = get_nitrate(datetime_object, lat, lon, pres[index], temp[index], psal[index], doxy[index])

    return output_suazade


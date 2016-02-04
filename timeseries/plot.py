# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import os.path as path
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

from datetime import datetime

from commons.helpers import is_number, get_date_string

class CoastEnum:
    coast, open_sea, everywhere = range(3)

    @staticmethod
    def valid(val):
        return val in range(3)

class SubBasinEnum:
    alb, sww, swe, nwm, tyr, adn, ads, aeg, ion, lev, med = range(11)

    @staticmethod
    def valid(val):
        return val in range(11)

class StatEnum:
    mean, std, p25, p50, p75 = range(5)

    @staticmethod
    def valid(val):
        return val in range(5)

def plot_from_files(file_list, varname, subbasin, coast=CoastEnum.open_sea, stat=StatEnum.mean, depth_index=0, fig=None, ax=None):
    """
    Plots a time series based on a list of file paths.

    Args:
        - *file_list*: a list of path strings for the files in the time series.
        - *varname*: name of the variable to plot.
        - *subbasin*: an element from SubBasinEnum.
        - *coast* (optional): an element from CoastEnum (default:
          CoastEnum.open_sea).
        - *stat* (optional): an element from StatEnum (default: StatEnum.mean).
        - *depth_index* (optional): the depth index z (default: 0).
        - *fig* (optional): an instance of matplotlib figure. A new one will be
          created if it is set to None (default: None).
        - *ax* (optional): an instance of matplotlib axes. A new one will be
          created if it is set to None (default: None).

    Returns: a matplotlib figure and axes object.
    """
    if (fig is None) or (ax is None):
        fig , ax = plt.subplots()
    plot_list = list()
    label_list = list()
    #For each file in file_list
    for f in file_list:
        #Get date string from file name
        _, ds = get_date_string(path.basename(f))
        #Create datetime object from date string
        dt = datetime.strptime(ds,'%Y%m%d')
        #Append the date to label_list
        label_list.append(dt)
        #Open it with netCDF4
        dset = netCDF4.Dataset(f)
        #Append the variable value to plot_list
        plot_list.append(dset[varname][subbasin, coast, depth_index, stat])
        #Close the file
        dset.close()
    #Plot data
    ax.plot(plot_list)
    #Set labels
    ax.set_xticklabels(label_list, rotation='vertical')
    return fig,ax

def plot_Hovmoeller_diagram(file_list, varname, subbasin, coast=CoastEnum.open_sea, stat=StatEnum.mean, depths=72, fig=None, ax=None):
    """
    Plots a time series Hovmoeller diagram.

    Args:
        - *file_list*: a list of path strings for the files in the time series.
        - *varname*: name of the variable to plot.
        - *subbasin*: an element from SubBasinEnum.
        - *coast* (optional): an element from CoastEnum (default:
          CoastEnum.open_sea).
        - *stat* (optional): an element from StatEnum (default: StatEnum.mean).
        - *depths* (optional): integer OR list OR Numpy array of depth values
          (default: 72).  If you pass a single integer it will be interpreted
          as the lenght of the depths array if you pass a list or Numpy array
          it will be used to set the labels and its lenght will be assumed as
          the lenght of the depths array.
        - *fig* (optional): an instance of matplotlib figure. A new one will be
          created if it is set to None (default: None).
        - *ax* (optional): an instance of matplotlib axes. A new one will be
          created if it is set to None (default: None).

    Returns: a matplotlib figure and axes object.
    """
    if (fig is None) or (ax is None):
        fig , ax = plt.subplots()
    dlabels = None
    if isinstance(depths, (int, long)):
        pass
    elif isinstance(depths, (list, tuple)):
        dlabels = np.array(depths)
        depths = len(dlabels)
    elif isinstance(depths, np.ndarray):
        #If 1D array
        if len(depths.shape) == 1:
            dlabels = np.copy(depths)
            depths = dlabels.shape[0]
        else:
            raise ValueError("Invalid depths argument")
    else:
        raise ValueError("Invalid depths argument")
    plotmat = np.zeros([depths, len(file_list)])
    label_list = list()
    #For each file
    for i,f in enumerate(file_list):
        #Get date string from file name
        _, ds = get_date_string(path.basename(f))
        #Create datetime object from date string
        dt = datetime.strptime(ds,'%Y%m%d')
        #Append the date to label_list
        label_list.append(dt)
        #Open it with netCDF4
        dset = netCDF4.Dataset(f)
        #Copy the data in the plot matrix
        plotmat[:,i] = dset[varname][subbasin, coast, 0:depths, stat]
        #Close the file
        dset.close()
    #Plot the matrix
    ax.imshow(plotmat)
    #Set labels
    ax.set_xticks(range(len(label_list)))
    ax.set_xticklabels(label_list, rotation='vertical')
    if not (dlabels is None):
        ticks = np.array(np.round(np.linspace(0, depths-1, num=8)), dtype=int)
        #ticks = [x * n for x in range(n)]
        labs = dlabels[ticks]
        ax.set_yticks(ticks)
        ax.set_yticklabels(labs)
    return fig,ax

if __name__ == "__main__":
    from glob import glob
    fl = sorted(glob('timeseries/*nc'))
    #fig,ax = plot_from_files(fl, 'O2o', SubBasinEnum.med)
    #plt.show()
    #fig,ax = plot_Hovmoeller_diagram(fl, 'O2o', SubBasinEnum.med)
    depths = [x * (5000 / 72) for x in range(72)]
    fig,ax = plot_Hovmoeller_diagram(fl, 'O2o', SubBasinEnum.med, depths=depths)
    plt.show()

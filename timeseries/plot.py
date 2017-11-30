# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import os.path as path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mpldates
import netCDF4
import pickle

from datetime import datetime

from commons.utils import is_number, get_date_string

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

def read_basic_info(stat_profile_file):
    '''
    Returns basin info
    '''
    ncIN = netCDF4.Dataset(stat_profile_file,'r')
    SUBLIST = str(ncIN.sub___list).split(", ")
    COASTLIST=str(ncIN.coast_list).split(", ")
    STAT_LIST=str(ncIN.stat__list).split(", ")
    ncIN.close()
    return SUBLIST, COASTLIST, STAT_LIST

def read_pickle_file(filename):
    print filename
    fid =open(filename,'r')
    [TIMESERIES,TL] = pickle.load(fid)
    fid.close()
    return TIMESERIES, TL


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
    ax.plot(label_list, plot_list)
    return fig,ax

def Hovmoeller_matrix(TIMESERIES, TL, depths, isub, icoast=0, istat=0):
    '''
    Arguments:
    TIMESERIES : numpy 4D array [nFrames,nSub,nCoast,jpk, nStat]
    TL         : Timelist object
    depths     : numpy array of nav_lev
    isub       : integer
    icoast     : integer
    istat      : integer
    '''
    ndepths=len(depths)
    plotmat = TIMESERIES[:,isub,icoast,0:ndepths,istat]
    nans_in_profile = np.isnan(plotmat[0,:])
    if not nans_in_profile.all(): # case of no nans
        first_nan = ndepths+1
    else:
        first_nan = np.argmax(nans_in_profile)

    plotmat = plotmat[:,:first_nan]
    depths  =  depths[  :first_nan]

    xlabel_list = mpldates.date2num(TL.Timelist)
    xs,ys = np.meshgrid(xlabel_list, depths)
    return plotmat.T, xs, ys

def Simple_timeseries(TIMESERIES,TL, iLev=0, iSub=0, iCoast=0, iStat=0):
    '''
    Retrives data for (x,y) plot
    '''
    return TL.Timelist, TIMESERIES[:,iSub,iCoast,iLev,iStat]


def Hovmoeller_matrix_from_SP(datetime_list, file_list, varname, subbasin, coast=CoastEnum.open_sea, stat=StatEnum.mean, depths=72):
    '''
    Reads directly from STAT_PROFILES
    '''
    dlabels=None
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
    #For each file
    for i,f in enumerate(file_list):
        #Open it with netCDF4
        dset = netCDF4.Dataset(f)
        #Copy the data in the plot matrix
        plotmat[:,i] = dset[varname][subbasin, coast, 0:depths, stat]
        #Close the file
        dset.close()
    xlabel_list = mpldates.date2num(datetime_list)
    xs,ys = np.meshgrid(xlabel_list, dlabels)
    return plotmat, xs, ys

def Hovmoeller_diagram(plotmat, xs,ys, fig=None, ax=None, shading='flat',vmin=None, vmax=None, cmap=None):
    if (fig is None) or (ax is None):
        fig , ax = plt.subplots()
    quadmesh = ax.pcolormesh(xs, ys, plotmat,shading=shading,vmin=vmin,vmax=vmax,cmap=cmap)
    #Inform matplotlib that the x axis is made by dates
    ax.xaxis_date()
    ax.invert_yaxis()
    return fig, ax, quadmesh
    
    


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

    Returns: a matplotlib Figure, Axes and QuadMesh object tuple.
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
    xlabel_list = list()
    #For each file
    for i,f in enumerate(file_list):
        #Get date string from file name
        _, ds = get_date_string(path.basename(f))
        #Create datetime object from date string
        dt = datetime.strptime(ds,'%Y%m%d')
        #Append the date to xlabel_list
        xlabel_list.append(dt)
        #Open it with netCDF4
        dset = netCDF4.Dataset(f)
        #Copy the data in the plot matrix
        plotmat[:,i] = dset[varname][subbasin, coast, 0:depths, stat]
        #Close the file
        dset.close()
    #Create the meshgrid
    xlabel_list = mpldates.date2num(xlabel_list)
    xs,ys = np.meshgrid(xlabel_list, dlabels)
    #Plot the matrix
    quadmesh = ax.pcolormesh(xs, ys, plotmat,shading='flat')
    #Inform matplotlib that the x axis is made by dates
    ax.xaxis_date()
    return fig, ax, quadmesh

if __name__ == "__main__":
    INPUTDIR="/Users/gbolzon/Documents/workspace/chain/postproc/"
    var="P_l"
    iSub=0
    iCoast=2
    iLev=0
    iStat=0
    plt.close('all')
    filename = INPUTDIR + var + ".pkl"
    TIMESERIES, TL = read_pickle_file(filename)

    fig,ax = plt.subplots()
    t,y = Simple_timeseries(TIMESERIES, TL, iLev, iSub, iCoast, iStat)
    ax.plot(t,y)
    fig.show()

    depths = np.arange(120)
    M, xs, ys = Hovmoeller_matrix(TIMESERIES,TL, depths, iSub, iCoast, iStat)
    fig, ax, quadmesh = Hovmoeller_diagram(M, xs, ys)

    fig.show()

    from glob import glob
    from commons.mask import Mask
    m = Mask('./layer_integral/meshmask.nc')
    fl = sorted(glob('timeseries/*nc'))
    #fig,ax = plot_from_files(fl, 'O2o', SubBasinEnum.med)
    #plt.show()
    #fig,ax = plot_Hovmoeller_diagram(fl, 'O2o', SubBasinEnum.med)
    depths = m.zlevels[0:30]
    fig,ax,im = plot_Hovmoeller_diagram(fl, 'O2o', SubBasinEnum.med, depths=depths)
    ax.invert_yaxis()
    fig.suptitle('O2o')
    fig.autofmt_xdate()
    plt.colorbar(im)
    plt.show()

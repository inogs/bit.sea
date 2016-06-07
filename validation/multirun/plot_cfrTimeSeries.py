import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Extracts statistcs from STAT_PROFILES
    and produces time series plots for run comparison

    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required =True,
                                help = ''' Output dir'''
                                )

    parser.add_argument(   '--indir', '-i',
                                type = str,
                                required =True,
                                help = ''' Input dir'''
                                )

    parser.add_argument(   '--satdir', '-s',
                                type = str,
                                required =False,
                                help = ''' Sat dir'''
                                )

    parser.add_argument(   '--pathtxt', '-t',
                                type = str,
                                required =True,
                                help = ''' Text of path of input files'''
                                )

    parser.add_argument(   '--cfrtype', '-c',
                                type = str,
                                required =True,
                                choices = ['nocoastda','daeffect','errmod','nutupt','vhany'],
                                help = ''' Type of comparison: nocoastda - runs without DAcoast (CR, 01, 06), daeffect - control run, run DAopensea and run DAcoast (CR, 01, 02), errmod - run DAopensea, run DAcoast and run DAcoast with modified model error (01, 02, 03), vhany - runs with and without vh anysotropic'''
                                )

    return parser.parse_args()


args = argument()

import os
import os.path as path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mpldates
import netCDF4

from glob import glob
from datetime import datetime

from commons.mask import Mask
from commons.utils import is_number, get_date_string


#VAR_LIST = ['P_n','P_p']
#VAR_LIST = ['ppn']
#VAR_LIST = ['N1p','N3n']
#VAR_LIST = ['P_i','P_c']
VAR_LIST = ['N1p']


DEPTH_LIST = [0, 50, 150]
# Non so come passarli

if args.cfrtype=='nocoastda':
   RUN_LIST = ['CR_COAST','DA_COAST_01','DA_COAST_06']
   LCOLOR_DICT={
               'CR_COAST': 'b',
               'DA_COAST_01': 'g',
               'DA_COAST_06': 'purple'
               }
   LEG_DICT={
            'CR_COAST': 'CR',
            'DA_COAST_01': '01',
            'DA_COAST_06': '06'
             }
if args.cfrtype=='daeffect':
   RUN_LIST = ['CR_COAST','DA_COAST_01','DA_COAST_02']
   LCOLOR_DICT={
               'CR_COAST': 'b',
               'DA_COAST_01': 'g',
               'DA_COAST_02': 'm'
               }
   LEG_DICT={
            'CR_COAST': 'CR',
            'DA_COAST_01': '01',
            'DA_COAST_02': '02'
            }
if args.cfrtype=='errmod':
   RUN_LIST = ['DA_COAST_01','DA_COAST_02','DA_COAST_03']
   LCOLOR_DICT={
               'DA_COAST_01': 'g',
               'DA_COAST_02': 'm',
               'DA_COAST_03': 'c'
               }
   LEG_DICT={
            'DA_COAST_01': '01',
            'DA_COAST_02': '02',
            'DA_COAST_03': '03'
            }
if args.cfrtype=='nutupt':
   RUN_LIST = ['DA_COAST_05','DA_COAST_03']
   LCOLOR_DICT={
               'DA_COAST_05': 'y',
               'DA_COAST_03': 'c'
               }
   LEG_DICT={
            'DA_COAST_05': '05',
            'DA_COAST_03': '03'
            }
if args.cfrtype=='vhany':
   RUN_LIST = ['DA_COAST_03','DA_COAST_04']
   LCOLOR_DICT={
               'DA_COAST_03': 'c',
               'DA_COAST_04': 'grey'
               }
   LEG_DICT={
            'DA_COAST_03': '03',
            'DA_COAST_04': '04'
            }



UNITS_DICT={
         'ppn' : 'mgC/m^3/d',
         'N1p' : 'mmol/m^3',
         'N3n' : 'mmol/m^3',
         'P_i' :'mgChl/m^3',
         'P_c' :'mgC/m^3',
         'P_n' :'mmolN/m^3',
         'P_p' :'mmolP/m^3'
         }

LSTYLE_DICT={
            0: '-',
           50: ':',
          150: '--' 
          }



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


SUB_LIST = ['alb', 'sww', 'swe', 'nwm', 'tyr', 'adn', 'ads', 'aeg', 'ion', 'lev', 'med']
COAST_LIST = ['coast', 'open_sea', 'everywhere']

maskfile = os.getenv('MASKFILE')
TheMask = Mask(maskfile)

def plot_from_runs(run_list, varname, subbasin, coast=CoastEnum.open_sea, stat=StatEnum.mean, fig=None, ax=None, satdir=None):
    """
    Plots a time series based on a list of file paths.

    Args:
        - *run_list*: a list runs to be compared
        - *varname*: name of the variable to plot.
        - *subbasin*: an element from SubBasinEnum.
        - *coast* (optional): an element from CoastEnum (default:
          CoastEnum.open_sea).
        - *stat* (optional): an element from StatEnum (default: StatEnum.mean).
        - *fig* (optional): an instance of matplotlib figure. A new one will be
          created if it is set to None (default: None).
        - *ax* (optional): an instance of matplotlib axes. A new one will be
          created if it is set to None (default: None).
        - *satdir* (optional): path of satellite directory for chl comparison

    Returns: a matplotlib figure and axes object.
    """
    if (fig is None) or (ax is None):
        fig , ax = plt.subplots()
    for dep in DEPTH_LIST:
        for run in RUN_LIST:
            run_path = args.indir + run + args.pathtxt
            file_list = sorted(glob(run_path + '/*nc'))
            plot_list = list()
            label_list = list()
            depth_index = TheMask.getDepthIndex(dep)
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
            ax.plot(label_list, plot_list, \
                    color=LCOLOR_DICT[run], linestyle=LSTYLE_DICT[dep], \
                    label=LEG_DICT[run]+ ' ' +str(dep)+'m')
        ax.set_ylabel(varname + ' [' + UNITS_DICT[varname] + ']')
        ax.legend(loc="best", \
                  fontsize='small', \
                  labelspacing=0, handletextpad=0,borderpad=0.1)
        ax.set_title(SUB_LIST[subbasin] + ' ' + \
                 COAST_LIST[coast]) 
    fig.autofmt_xdate()
    return fig,ax


inpsat=None
for isub in range(len(SUB_LIST)):
    print(SUB_LIST[isub])
    for var in VAR_LIST:
        if var=='P_i': inpsat=args.satdir 
        print(var)
        fig1 = plt.figure(num=None,facecolor='w', \
                        edgecolor='k',figsize=[6.5,8.5]) 
        ax1  = fig1.add_subplot(2,1,1)
        plot_from_runs(RUN_LIST, var, isub, \
                   coast = CoastEnum.open_sea, \
                   fig=fig1,ax=ax1,satdir=inpsat)
        ax2  = fig1.add_subplot(2,1,2)
        plot_from_runs(RUN_LIST, var, isub, \
                   coast = CoastEnum.coast, \
                   fig=fig1,ax=ax2,satdir=inpsat)

        outfile = args.outdir + var + '_' + \
                  args.cfrtype + '_' + \
                  SUB_LIST[isub] + '.png'
        plt.savefig(outfile)
      #  plt.show()
        plt.close(fig1)

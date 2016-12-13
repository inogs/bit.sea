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

    parser.add_argument(   '--cfrtype', '-c',
                                type = str,
                                required =True,
                                choices = ['nocoastda','daeffect','errmod','nutupt','vhany','theta','alpha','errsat','short','corrl','dacoast'],
                                help = ''' Type of comparison: nocoastda - runs without DAcoast (CR, 01, 06), daeffect - control run, run DAopensea and run DAcoast (CR, 01, 02), errmod - run DAopensea, run DAcoast and run DAcoast with modified model error (01, 02, 03), nutupt - CR, run without and with modified formulation of nutrient uptake (CR, 01, 06), vhany - runs with and without vh anysotropic, theta - different chl:c imposed in DA (01, 07, 08), alpha - effect of alpha +5% -10% (07, 09, 10), errsat - augmented errsat (CR, 08, 11, 12), short - short tests on errsat, corrl - correlation radius 10km (12, 13), dacoast - effect of coastal DA (CR, 13, 14)'''
                                )

    return parser.parse_args()


args = argument()

import os
import os.path as path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mpldates
import pickle

from glob import glob



#VAR_LIST = ['P_n','P_p']
#VAR_LIST = ['ppn']
#VAR_LIST = ['N1p']
#VAR_LIST = ['N3n']
VAR_LIST = ['P_i','P_c','P_n','P_p']
VAR_LIST = ['ppn','N1p','N3n']


DEPTH_LIST = [0, 50, 150]
# Non so come passarli

if args.cfrtype=='nocoastda':
   RUN_LIST = ['CR_COAST','DA_COAST_01','DA_COAST_06','DA_COAST_07']
   LCOLOR_DICT={
               'CR_COAST': 'b',
               'DA_COAST_01': 'g',
               'DA_COAST_06': 'purple',
               'DA_COAST_07': 'y'
               }
   LEG_DICT={
            'CR_COAST': 'CR',
            'DA_COAST_01': '01',
            'DA_COAST_06': '06',
            'DA_COAST_07': '07'
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
   RUN_LIST = ['CR_COAST','DA_COAST_01','DA_COAST_06']
   LCOLOR_DICT={
               'CR_COAST'   : 'b',
               'DA_COAST_01': 'g',
               'DA_COAST_06': 'purple'
               }
   LEG_DICT={
            'CR_COAST'   : 'CR',
            'DA_COAST_01': '01',
            'DA_COAST_06': '06'
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
if args.cfrtype=='theta':
   RUN_LIST = ['CR_COAST','DA_COAST_06','DA_COAST_07','DA_COAST_08']
   LCOLOR_DICT={
               'CR_COAST': 'b',
               'DA_COAST_06': 'purple',
               'DA_COAST_07': 'y',
               'DA_COAST_08': 'grey'
               }
   LEG_DICT={
            'CR_COAST': 'CR',
            'DA_COAST_06': '06',
            'DA_COAST_07': '07',
            'DA_COAST_08': '08'
            }

if args.cfrtype=='alpha':
   RUN_LIST = ['DA_COAST_07','DA_COAST_09','DA_COAST_10']
   LCOLOR_DICT={
               'DA_COAST_07': 'y',
               'DA_COAST_09': 'g',
               'DA_COAST_10': 'b'
               }
   LEG_DICT={
            'DA_COAST_07': '07',
            'DA_COAST_09': '09',
            'DA_COAST_10': '10'
            }

if args.cfrtype=='errsat':
   RUN_LIST = ['CR_COAST','DA_COAST_08','DA_COAST_11','DA_COAST_12']
   LCOLOR_DICT={
               'CR_COAST': 'b',
               'DA_COAST_08': 'grey',
               'DA_COAST_11': 'r',
               'DA_COAST_12': 'c'
               }
   LEG_DICT={
            'CR_COAST': 'CR',
            'DA_COAST_08': '08',
            'DA_COAST_11': '11',
            'DA_COAST_12': '12'
            }

if args.cfrtype=='short':
   RUN_LIST = ['DA_COAST_08','DA_COAST_11','VARSAT50','VARSAT100']
   LCOLOR_DICT={
               'DA_COAST_08': 'grey',
               'DA_COAST_11': 'r',
               'VARSAT50': 'k',
               'VARSAT100': 'g'
               }
   LEG_DICT={
            'DA_COAST_08': '08',
            'DA_COAST_11': '11',
            'VARSAT50': '+50%',
            'VARSAT100': '+100%'
            }

if args.cfrtype=='corrl':
   RUN_LIST = ['DA_COAST_12','DA_COAST_13']
   LCOLOR_DICT={
               'DA_COAST_12': 'r',
               'DA_COAST_13': 'grey',
               'DA_COAST_14': 'g'
               }
   LEG_DICT={
            'DA_COAST_12': '12',
            'DA_COAST_13': '13',
            'DA_COAST_14': '14'
            }


if args.cfrtype=='dacoast':
   RUN_LIST = ['CR_COAST','DA_COAST_13','DA_COAST_14']
   LCOLOR_DICT={
               'CR_COAST'   : 'b',
               'DA_COAST_13': 'grey',
               'DA_COAST_14': 'g'
               }
   LEG_DICT={
            'CR_COAST'   : 'CR',
            'DA_COAST_13': '13',
            'DA_COAST_14': '14'
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
           50: '--',
          150: ':' 
          }



SUB_LIST = ['alb', 'sww', 'swe', 'nwm', 'tyr', 'adn', 'ads', 'aeg', 'ion', 'lev', 'med']


def plot_from_runs(run_list, varname, subbasin, coast='open', fig=None, ax=None, satdate=None, satchl=None):
    """
    Plots a time series based on a list of file paths.

    Args:
        - *run_list*: a list runs to be compared
        - *varname*: name of the variable to plot
        - *subbasin*: a subbasin nam
        - *coast* (optional): 'open' or 'coast'  (default:'open').
        - *fig* (optional): an instance of matplotlib figure. A new one will be
          created if it is set to None (default: None).
        - *ax* (optional): an instance of matplotlib axes. A new one will be
          created if it is set to None (default: None).
        - *satdate* (optional): dates of satellite observations
        - *satchl* (optional): values of satellite observations

    Returns: a matplotlib figure and axes object.
    """
    if (fig is None) or (ax is None):
        fig , ax = plt.subplots()
    for run in RUN_LIST:
        txtfile = varname + '_' + \
                  run + '_' + \
                  subbasin + \
                  coast + '.pkl'
        print(txtfile)
        file_pkl = args.indir + txtfile
        fid = open(file_pkl)
        LIST = pickle.load(fid)
        fid.close()
        label_list = LIST[0]
        for idep,dep in enumerate(DEPTH_LIST):
            plot_list = LIST[1+idep]
            if (dep==0):
                    ax.plot(label_list, plot_list, \
                    color=LCOLOR_DICT[run], linestyle=LSTYLE_DICT[dep], \
                    label=LEG_DICT[run]+ ' ' +str(dep)+'m')
            else:
                    ax.plot(label_list, plot_list, \
                    color=LCOLOR_DICT[run], linestyle=LSTYLE_DICT[dep])
    if satdate:
       ax.plot(satdate,satchl,'xk', label='SAT')

    ax.set_ylabel(varname + ' [' + UNITS_DICT[varname] + ']')
    ax.legend(loc="best", \
              fontsize='small', \
              labelspacing=0, handletextpad=0,borderpad=0.1)
    ax.set_title(subbasin + ' ' + coast) 
    fig.autofmt_xdate()

    return fig,ax


MY_YEAR ='2013'
print('Year for sat hard coded (only for P_i): ' + MY_YEAR + '!')

datatypes = ['open','coast']
DICsat_file = {
         'coast': args.satdir + 'satstats.coast' + MY_YEAR + '.pkl',
         'open' : args.satdir + 'satstats.open' + MY_YEAR + '.pkl'
              }

LISTsat = [i for i in range(3)]
for it,typed in enumerate(datatypes):
    fid = open(DICsat_file[typed])
    LIST = pickle.load(fid)
    fid.close()
    LISTsat[0] = LIST[0]
    LISTsat[it+1] = LIST[1]


for isub,sub in enumerate(SUB_LIST):
    print(sub)
    for var in VAR_LIST:
        print(var)
        if (var is 'P_i'):
           datesat=LISTsat[0]
           chlsat=LISTsat[1]
           chlsatsub =chlsat[isub,:]
        else: 
           datesat = None
           chlsatsub = None
        fig1 = plt.figure(num=None,facecolor='w', \
                        edgecolor='k',figsize=[6.5,8.5]) 
        ax1  = fig1.add_subplot(2,1,1)
        plot_from_runs(RUN_LIST, var, sub, \
                   coast = 'open', \
                   fig=fig1,ax=ax1,satdate=datesat,satchl=chlsatsub)
        if (var is 'P_i'):
           chlsat=LISTsat[1]
        else: 
           chlsat = None
        ax2  = fig1.add_subplot(2,1,2)
        plot_from_runs(RUN_LIST, var, sub, \
                   coast = 'coast', \
                   fig=fig1,ax=ax2,satdate=datesat,satchl=chlsatsub)

        outfile = args.outdir + var + '_' + \
                  args.cfrtype + '_' + \
                  sub + '.png'
        plt.savefig(outfile)
      #  plt.show()
        plt.close(fig1)

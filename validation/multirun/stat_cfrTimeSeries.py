import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Produces icomparison between time series

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

    parser.add_argument(   '--cfrtype', '-c',
                                type = str,
                                required =True,
                                choices = ['nocoastda','daeffect','errmod','nutupt','vhany','theta','alpha','errsat','short','corrl','dacoast'],
                                help = ''' Type of comparison: nocoastda - runs without DAcoast (CR, 01, 06), daeffect - control run, run DAopensea and run DAcoast (CR, 01, 02), errmod - run DAopensea, run DAcoast and run DAcoast with modified model error (01, 02, 03), nutupt - CR, run without and with modified formulation of nutrient uptake (CR, 01, 06), vhany - runs with and without vh anysotropic, theta - different chl:c imposed in DA (01, 07, 08), alpha - effect of alpha +5% -10% (07, 09, 10), errsat - augmented errsat (CR, 08, 11, 12), short - short tests on errsat, corrl - correlation radius 10km (08, 13), dacoast - effect of coastal DA (CR, 13, 14)'''
                                )

    return parser.parse_args()


args = argument()

import os
import os.path as path
import numpy as np
import pickle

from glob import glob



VAR_LIST = ['ppn','N1p','N3n','P_i','P_c','P_n','P_p']

DEPTH_LIST = [0, 50, 150]

SUB_LIST = ['alb', 'sww', 'swe', 'nwm', 'tyr', 'adn', 'ads', 'aeg', 'ion', 'lev', 'med']
 
selectvar = 'P_i'
selectsub = 'med'


DiffStat = {}

for sub in SUB_LIST:
    DiffStat[sub] = {}
    for var in VAR_LIST:
        DiffStat[sub][var] = {}

if args.cfrtype=='nocoastda':
   RUN_LIST = ['CR_COAST','DA_COAST_01','DA_COAST_06','DA_COAST_07']

if args.cfrtype=='daeffect':
   RUN_LIST = ['CR_COAST','DA_COAST_01','DA_COAST_02']

if args.cfrtype=='errmod':
   RUN_LIST = ['DA_COAST_01','DA_COAST_02','DA_COAST_03']

if args.cfrtype=='nutupt':
   RUN_LIST = ['CR_COAST','DA_COAST_01','DA_COAST_06']

if args.cfrtype=='vhany':
   RUN_LIST = ['DA_COAST_03','DA_COAST_04']

if args.cfrtype=='theta':
   RUN_LIST = ['DA_COAST_01','DA_COAST_07','DA_COAST_08']

if args.cfrtype=='alpha':
   RUN_LIST = ['DA_COAST_07','DA_COAST_09','DA_COAST_10']

if args.cfrtype=='errsat':
   RUN_LIST = ['CR_COAST','DA_COAST_08','DA_COAST_11','DA_COAST_12']

if args.cfrtype=='short':
   RUN_LIST = ['DA_COAST_08','DA_COAST_11','VARSAT50','VARSAT100']

if args.cfrtype=='corrl':
   RUN_LIST = ['DA_COAST_12','DA_COAST_13']

if args.cfrtype=='dacoast':
   RUN_LIST = ['CR_COAST','DA_COAST_13','DA_COAST_14']


UNITS_DICT={
         'ppn' : 'mgC/m^3/d',
         'N1p' : 'mmol/m^3',
         'N3n' : 'mmol/m^3',
         'P_i' :'mgChl/m^3',
         'P_c' :'mgC/m^3',
         'P_n' :'mmolN/m^3',
         'P_p' :'mmolP/m^3'
         }


def rmsd(list1,list2):
    RMSDvalue = []
    if list1:
       listarr1 = np.array(list1)
       listarr2 = np.array(list2)
       RMSDvalue = (np.mean((listarr2-listarr1)**2))**.5
    return RMSDvalue

def stat_from_runs(run_list, varname, subbasin, coast='open'):
    """
    Plots a time series based on a list of file paths.

    Args:
        - *run_list*: a list runs to be compared.
        - *varname*: name of the variable tobe compared. 
        - *subbasin*: a subbasin name.
        - *coast* (optional): 'open' or 'coast'  (default:'open').

    Returns: statistics of differences. 
    """
    rmsd_values = {}
    times_list = {}
    for dep in DEPTH_LIST:
        rmsd_values[dep] = []
        times_list[dep] = []
    for irun,run in enumerate(RUN_LIST):
        txtfile = varname + '_' + \
                  run + '_' + \
                  subbasin + \
                  coast + '.pkl'
        file_pkl = args.indir + txtfile
        fid = open(file_pkl)
        LIST = pickle.load(fid)
        fid.close()
        for idep,dep in enumerate(DEPTH_LIST):
            rmsd_values[dep].append(rmsd(times_list[dep],LIST[1+idep]))
            times_list[dep] = LIST[1+idep]
            

    return rmsd_values



for sub in SUB_LIST:
    print(sub)
    for var in VAR_LIST:
        print(var)
        DiffStat[sub][var] = stat_from_runs(RUN_LIST, var, sub, \
                   coast = 'open')
        DiffStat[sub][var] = stat_from_runs(RUN_LIST, var, sub, \
                   coast = 'coast')

outfile = args.outdir +  'statdiff_' + \
          args.cfrtype + '.pkl'

fid = open(outfile,'wb')
pickle.dump(DiffStat,fid)
fid.close()

print(DiffStat[selectsub][selectvar])
print(RUN_LIST)

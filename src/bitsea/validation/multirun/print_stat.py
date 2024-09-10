import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Produces icomparison between time series

    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--indir', '-i',
                                type = str,
                                required =True,
                                help = ''' Input dir'''
                                )

    parser.add_argument(   '--cfrtype', '-c',
                                type = str,
                                required =True,
                                choices = ['nocoastda','daeffect','errmod','nutupt','vhany','theta','alpha','errsat','short','corrl'],
                                help = ''' Type of comparison: nocoastda - runs without DAcoast (CR, 01, 06), daeffect - control run, run DAopensea and run DAcoast (CR, 01, 02), errmod - run DAopensea, run DAcoast and run DAcoast with modified model error (01, 02, 03), nutupt - CR, run without and with modified formulation of nutrient uptake (CR, 01, 06), vhany - runs with and without vh anysotropic, theta - different chl:c imposed in DA (01, 07, 08), alpha - effect of alpha +5% -10% (07, 09, 10), errsat - augmented errsat (CR, 01, 11, 12), short - short tests on errsat, corrl - correlation radius 10km (08, 13)'''
                                )

    return parser.parse_args()


args = argument()

import numpy as np
import pickle


selectvar = 'P_i'
selectsub = 'med'


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
   RUN_LIST = ['CR_COAST','DA_COAST_01','DA_COAST_11','DA_COAST_12']

if args.cfrtype=='short':
   RUN_LIST = ['DA_COAST_08','DA_COAST_11','VARSAT50','VARSAT100']

if args.cfrtype=='corrl':
   RUN_LIST = ['DA_COAST_11','DA_COAST_13']


UNITS_DICT={
         'ppn' : 'mgC/m^3/d',
         'N1p' : 'mmol/m^3',
         'N3n' : 'mmol/m^3',
         'P_i' :'mgChl/m^3',
         'P_c' :'mgC/m^3',
         'P_n' :'mmolN/m^3',
         'P_p' :'mmolP/m^3'
         }



pklfile = args.indir +  'statdiff_' + \
          args.cfrtype + '.pkl'

fid = open(pklfile)
Stat = pickle.load(fid)
fid.close()

print(selectsub,selectvar)
print(RUN_LIST)
print(Stat[selectsub][selectvar])

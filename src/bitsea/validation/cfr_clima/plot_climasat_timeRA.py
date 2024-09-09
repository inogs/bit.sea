import numpy as np
import argparse
import pickle
import matplotlib.pyplot as plt
import scipy.io.netcdf as NC
from datetime import datetime, timedelta
#from Sat import SatManager as Sat
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.dataextractor import DataExtractor
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2
#from postproc import masks
from commons.utils import addsep
from profilerplot import plotNAMESLIST
#import os


def argument():
    parser = argparse.ArgumentParser(description = '''
    Plot model and sat climatology
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' FIGURES  directory'''
                                )


    parser.add_argument(   '--insat', '-s',
                                type = str,
                                required = True,
                                help = ''' Satellite directory'''
                                )


    parser.add_argument(   '--timeperiod', '-t',
                                type = int,
                                required = False,
                                help = ''' Time period for the plot (optional):
                                           number of months.
                                           If not specified, considers all the time period'''
                                )

    return parser.parse_args()


args = argument()


INSAT   = addsep(args.insat)
OUTDIR  = addsep(args.outdir)
if args.timeperiod is not None:
    monthWindow = args.timeperiod



LIST = {}
for listname in plotNAMESLIST:
    filestats = INSAT + '/' + listname +'.pkl'
    fid = open(filestats,'r')
    LIST[listname] = pickle.load(fid)
    fid.close()



# Plots 
plt.close('all')

areaSubLIST = {}
areaSubLIST['West1'] = ['alb','swm1','swm2']
areaSubLIST['West2'] = ['nwm','tyr1','tyr2']
areaSubLIST['Adr']   = ['adr1','adr2']
areaSubLIST['Ion']   = ['ion1','ion2','ion3']
areaSubLIST['Lev']   = ['lev1','lev2','lev3','lev4']
areaSubLIST['MED']   = ['med']

DICTsubcolor = {
    'alb': 'orange',
    'swm1': 'purple',
    'swm2': 'green',
    'nwm': 'cyan',
    'tyr1': 'violet',
    'tyr2': 'orange',
    'adr1': 'brown',
    'adr2': 'green',
    'ion1': 'darkolivegreen',
    'ion2': 'darkorange',
    'ion3': 'blue',
    'lev1': 'gold',
    'lev2': 'pink',
    'lev3': 'olive',
    'lev4': 'magenta',
    'med' : 'cornflowerblue',
}

indMask = {
    'Everywhere': 0,
    'Opensea': 1,
}
plt.close('all')


for txtarea in ['Everywhere','Opensea']:

    indexm = indMask[txtarea]
    for area in areaSubLIST.keys():
        print area + ' --- ' + txtarea
        fig,axs = plt.subplots(len(areaSubLIST[area]),1,sharex=True,figsize=(12,8))
        indax = 0
        for sub in V2.P:
            if sub.name in areaSubLIST[area]:
                if len(areaSubLIST[area])>1:
                    plt.sca(axs[indax])
                    indax += 1
                else:
                    plt.sca(axs)
                for listname in plotNAMESLIST:
                    plotDates = LIST[listname][0]
                    climaMean = LIST[listname][2]
                    satMean = LIST[listname][1]
                    climaStd = LIST[listname][3]
                    plt.vlines(plotDates, \
                        climaMean[sub.name][indexm,:]-2*climaStd[sub.name][indexm,:], \
                        climaMean[sub.name][indexm,:]+2*climaStd[sub.name][indexm,:], \
                        color='lightgrey')
                    plt.plot(plotDates,climaMean[sub.name][indexm,:],'-', \
                        color=DICTsubcolor[sub.name], \
                        lw=5,alpha=0.5, \
                        label='Satellite climatology ', \
                        )
                    plt.plot(plotDates,satMean[sub.name][indexm,:],'.', \
                        color='grey', \
                        label='Satellite assimilated', \
                        )
                    plt.ylabel('mg chl/m^3')
                    plt.grid()
                    plt.title(txtarea + ' ' + sub.name)
                    plt.legend(loc='best')
                    if args.timeperiod is not None:
                        _,xend = plt.xlim()
                        xstart = xend - timedelta(days=monthWindow*30)
                        plt.xlim(xstart,xend)


            plt.savefig(OUTDIR + 'satseries_clima' + area + '_' + txtarea + '.png')





plt.show(block=False)



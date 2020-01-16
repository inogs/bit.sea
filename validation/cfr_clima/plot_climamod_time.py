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

    parser.add_argument(   '--satclimdir', '-c',
                                type = str,
                                required = True,
                                help = ''' Climatology statistics directory'''
                                )

    parser.add_argument(   '--insat', '-s',
                                type = str,
                                required = True,
                                help = ''' Satellite directory'''
                                )

    parser.add_argument(   '--inmod', '-i',
                                type = str,
                                required = True,
                                help = ''' Model chlorophyll output directory'''
                                )

    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = ''' File of the model mask
                                           '''
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


CLIMDIR = addsep(args.satclimdir)
INDIR   = addsep(args.inmod)
INSAT   = addsep(args.insat)
OUTDIR  = addsep(args.outdir)
if args.timeperiod is not None:
    monthWindow = args.timeperiod

TheMask = Mask(args.maskfile)
mask200 = TheMask.mask_at_level(200)

print 'Subbasin masks'
submask = {}
subm_open = {}
subnames = ''
for sub in V2.Pred:
    print sub.name
    submask[sub.name] = SubMask(sub,maskobject=TheMask).mask[0,:,:]
    subm_open[sub.name] = submask[sub.name] | (mask200==True)
    subnames += sub.name + ', '

subnames += 'med, '
submask['med'] = np.zeros_like(submask[sub.name])
for sub in V2.Pred:
    submask['med'] = submask['med'] | submask[sub.name]
subm_open['med'] = submask['med'] | (mask200==True)



Timestart = "19000101"
#Timestart = "20160101"
Time__end = "20500101"

TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TLmodel = TimeList.fromfilenames(TI, INDIR, \
            "chl.*nc", prefix='chl.')

if args.timeperiod is not None:
    datetEnd   = TLmodel.Timelist[-1]
    datetStart = datetEnd - timedelta(days=monthWindow*30)
    Timestart = datetime.strftime(datetStart,"%Y%m%d")
    Time__end = datetime.strftime(datetEnd,"%Y%m%d")
    TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
    TLmodel = TimeList.fromfilenames(TI, INDIR, \
            "chl.*nc", prefix='chl.')

Ndates = TLmodel.nTimes


# Check on sublist of clima and sat statistics
filestats0 = CLIMDIR + 'statsday001_clim.nc'
NCclim = NC.netcdf_file(filestats0,'r')
# filesat0 = TLsat.filelist[0]
# NCsat = NC.netcdf_file(filestats0,'r')
if (not(NCclim.sublist==subnames)): # or (not(NCsat.sublist==subnames)):
    print '''List of subbasin used for climatology statistics 
             not equal to the actual one'''
    print 'CLIMA list:  ' + NCclim.sublist
    # print 'SAT list:    ' + NCsat.sublist
    print 'ACTUAL list: ' + subnames
    print ' ---  EXIT!  ---'
    NCclim.close()
    import sys
    sys.exit()

NCclim.close()


plotDates = []
modelMean = {}
#modelStd = {}
climaMean = {}
climaStd = {}
for sub in V2.P:
    modelMean[sub.name] = np.zeros((2,Ndates))
    #modelStd[sub.name] = np.zeros(Ndates)
    climaMean[sub.name] = np.zeros((2,Ndates))
    climaStd[sub.name] = np.zeros((2,Ndates))

for ifile,filemod in enumerate(TLmodel.filelist):
    datefile = TLmodel.Timelist[ifile]
    plotDates.append(datefile)
    print ' Reading model file ' + np.str(ifile+1) + ' of ' + np.str(Ndates)
    print '   ... date ' + datetime.strftime(datefile,'%Y-%m-%d')
    De = DataExtractor(TheMask,filename=filemod,varname='lchlm',dimvar=2)
    surfchl = De.filled_values
    
    julianday = int(datefile.strftime('%j'))
    if julianday==366:
        julianday = 365
    fileclima = CLIMDIR + 'statsday' + '%03d' %julianday + '_climAll.nc'
    NCclim = NC.netcdf_file(fileclima,'r')
    climastats = NCclim.variables['SUBstatistics'].data.copy()
    NCclim.close()

    fileclima = CLIMDIR + 'statsday' + '%03d' %julianday + '_climOpen.nc'
    NCclim = NC.netcdf_file(fileclima,'r')
    climast_open = NCclim.variables['SUBstatistics'].data.copy()
    NCclim.close()

    for isub,sub in enumerate(V2.P):
        modelMean[sub.name][0,ifile] = np.nanmean(surfchl[submask[sub.name]])
        modelMean[sub.name][1,ifile] = np.nanmean(surfchl[subm_open[sub.name]])
        #modelStd[sub.name][ifile] = np.nanstd(surfchl[submask[sub.name]])
        climaMean[sub.name][0,ifile] = climastats[isub,0]
        climaMean[sub.name][1,ifile] = climast_open[isub,0]
        climaStd[sub.name][0,ifile] = climastats[isub,9] #average std in the subbasin  
        climaStd[sub.name][1,ifile] = climast_open[isub,9] #average std in the subbasin  



TLsat = TimeList.fromfilenames(TI, INSAT, \
            "*.nc",prefix='',dateformat='%Y%m%d')
Ndsat = TLsat.nTimes

satMean = {}
datepsat = []
for sub in V2.P:
    satMean[sub.name] = np.zeros((2,Ndsat))

for ifile,filesat in enumerate(TLsat.filelist):
    datefile = TLsat.Timelist[ifile]
    datepsat.append(datefile)
    print ' Reading satellite file ' + np.str(ifile+1) + ' of ' + np.str(Ndsat)
    print '   ... date ' + datetime.strftime(datefile,'%Y-%m-%d')
    De = DataExtractor(TheMask,filename=filesat,varname='lchlm',dimvar=2)
    satchl = De.filled_values
    for isub,sub in enumerate(V2.P):
        satMean[sub.name][0,ifile] = np.nanmean(satchl[submask[sub.name]])
        satMean[sub.name][1,ifile] = np.nanmean(satchl[subm_open[sub.name]])


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
def plotsclima(climaMean,climaStd,satMean,modelMean,txtarea):
    plt.close('all')
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
                plt.vlines(plotDates, \
                        climaMean[sub.name][indexm,:]-2*climaStd[sub.name][indexm,:], \
                        climaMean[sub.name][indexm,:]+2*climaStd[sub.name][indexm,:], \
                        color='lightgrey')
                plt.plot(plotDates,climaMean[sub.name][indexm,:],'-', \
                    color=DICTsubcolor[sub.name], \
                    lw=5,alpha=0.5, \
                    label='Satellite climatology ', \
                    )
                # plt.plot(plotDates,climaMean[sub.name][indexm-1,:],'-', \
                #     color=DICTsubcolor[sub.name], \
                #     # lw=5,alpha=0.5, \
                #     label='Satellite climatology ', \
                #     )
                plt.plot(datepsat,satMean[sub.name][indexm,:],'.', \
                    color='grey', \
                    label='Satellite assimilated', \
                    )
                plt.plot(plotDates,modelMean[sub.name][indexm,:],'.',
                    color=DICTsubcolor[sub.name], \
                    label='Model ', \
                    )
                plt.ylabel('mg chl/m^3')
                plt.grid()
                plt.title(txtarea + ' ' + sub.name)
                plt.legend(loc='best')
        
        plt.savefig(OUTDIR + 'model_vs_clima' + area + '_' + txtarea + '.png')


plotsclima(climaMean,climaStd,satMean,modelMean,'Everywhere')
plotsclima(climaMean,climaStd,satMean,modelMean,'Opensea')



plt.show(block=False)



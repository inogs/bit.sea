import numpy as np
import argparse
import pickle
import matplotlib.pyplot as plt
#from Sat import SatManager as Sat
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins import V2
#from postproc import masks
from commons.utils import addsep
#import os


def argument():
    parser = argparse.ArgumentParser(description = '''
    Plot results of blacklisting comparing sat values with climatology
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' FIGURES  directory'''
                                )

    parser.add_argument(   '--indir', '-i',
                                type = str,
                                required = True,
                                help = ''' PKL  directory'''
                                )
    
    parser.add_argument(   '--submaskdir', '-m',
                                type = str,
                                required = True,
                                help = ''' Directory of submasks for satellite 
                                            grid (output of sat_indsub.py)'''
                                )

    return parser.parse_args()


args = argument()


INDIR   = addsep(args.indir)
OUTDIR  = addsep(args.outdir)
SUBMASKDIR  = addsep(args.submaskdir)

Timestart="19501231"
Time__end="20500101"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")


plt.close('all')

#plt.figure(num=1)

for sub in V2.P:
#for sub in [V2.med]: #V2.P:
    filemasksub = SUBMASKDIR + 'masksub.' + sub.name + '.npy'
    masksub = np.load(filemasksub)
    Ntotsub = np.sum(masksub)

ss1 = {}
for masktype in ['ORIG','CHECK','ALLP']:
    print(masktype)
    TLtime = TimeList.fromfilenames(TI, INDIR , \
                "stats_time*" + masktype + "*.pkl", \
                prefix='stats_time',dateformat='%Y%m%d')
    TLclim = TimeList.fromfilenames(TI, INDIR , \
                "stats_clim*" + masktype + "*.pkl", \
                prefix='stats_clim',dateformat='%Y%m%d')
    datesplot = TLtime.Timelist

    statsSUB_time = {}
    statsSUB_clim = {}
    for sub in V2.P: #[V2.med]:
        statsSUB_time[sub.name] = np.zeros((TLtime.nTimes,8))
        statsSUB_clim[sub.name] = np.zeros((TLclim.nTimes,12))

    for ifile,timefile in enumerate(TLtime.filelist):
        datesplot[ifile].replace(hour=23)
        datefile = TLtime.Timelist[ifile]
        #print(timefile)
        fid = open(timefile,'r')
        stats_timeday = pickle.load(fid)
        fid.close()

        iclim = TLclim.find(datefile)
        climfile = TLclim.filelist[iclim]
        fid = open(climfile,'r')
        stats_climday = pickle.load(fid)
        fid.close()

        
        for sub in V2.P: #[V2.med]: 
            statsSUB_time[sub.name][ifile,:] = stats_timeday[sub.name]
            statsSUB_clim[sub.name][ifile,:] = stats_climday[sub.name]

    for sub in V2.P:
        print(sub)
    #for sub in [V2.med]: #V2.P:
        fig,axs = plt.subplots(4,1,sharex=True,figsize=[8,11])

        nax = 0
        plt.sca(axs[0])
        plt.title('Cholophyll statistics over satellite coverage [ngchl/m^3]')
        plt.suptitle(sub.name.upper() + ' ' + masktype)

        plt.vlines(TLclim.Timelist, \
                    0.01,
                    statsSUB_clim[sub.name][:,3], \
                    color='grey')
        
        plt.plot(TLclim.Timelist, \
                    statsSUB_clim[sub.name][:,0],'o', \
                    color='grey',label='Mean climatology')

        ss1[masktype] = statsSUB_time[sub.name][:,0]
        plt.plot(datesplot, \
                    statsSUB_time[sub.name][:,0],'o', \
                    color='y',label='Mean obs')

        plt.legend(loc='best')
        plt.ylim(0.0,0.6)
        plt.grid()


        plt.sca(axs[1])
        plt.title('Mean |obs-meanclim|/stdclim')

        plt.plot(datesplot,statsSUB_time[sub.name][:,3])

        # axs[1].set_yscale('log')
        plt.ylim(0,5)
        plt.grid()


        plt.sca(axs[2])
        plt.title('Excluded observations [grid points/sub grid points with obs]')

        plt.bar(datesplot,(statsSUB_time[sub.name][:,6]+statsSUB_time[sub.name][:,7])/ \
                           statsSUB_time[sub.name][:,5])

        plt.grid()



        plt.sca(axs[3])
        plt.title('Coverage [grid points with obs / total sub grid points] ' + \
                   np.str(Ntotsub))

        plt.bar(datesplot,statsSUB_time[sub.name][:,5]/Ntotsub)

        plt.grid()

        plt.savefig('testfig' + masktype + '_' + sub.name + '.png')


    #plt.show(block=False)





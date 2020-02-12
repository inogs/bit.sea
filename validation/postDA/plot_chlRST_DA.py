import numpy as np
from commons.submask import SubMask
from basins import V2 as OGS
import pickle
import datetime
import matplotlib.pyplot as plt


RUNlist = ['RUN_TEST','TRANSITION']
RUNlist = ['RUN_2017Dard','TRANSITION']
masktype = 'open'
masktype = 'everywhere'

SATdir = '/gpfs/scratch/userexternal/ateruzzi/ELAB_CFR_24/SAT/'
PKLdir = '/gpfs/scratch/userexternal/ateruzzi/ELAB_CFR_24/RST_DA/'
OUTdir = '/gpfs/scratch/userexternal/ateruzzi/ELAB_CFR_24/FIGURES_DA/'

Nsub = len(OGS.P.basin_list)

times_b = {}
values_b = {}
times_a = {}
values_a = {}
tt = {}
indt = {}

fid = open(SATdir + 'meansat.pkl','r')
satLIST = pickle.load(fid)
fid.close()
sattimes = satLIST[0]
satvalues = satLIST[1][masktype]
satcov = satLIST[2][masktype]

fid = open(SATdir + 'pointsubbasin.pkl','r')
pointsub = pickle.load(fid)
npoints = pointsub[masktype]
fid.close()



misftimes = {}
misfvalues = {}
RMSvalues = {}
for run in RUNlist:
    fid = open(PKLdir + '/../MISF_DA/' + run + '/meanmisf.pkl','r')
    misfLIST = pickle.load(fid)
    fid.close()
    misftimes[run] = misfLIST[0]
    misfvalues[run] = misfLIST[1][masktype]
    RMSvalues[run] = misfLIST[2][masktype]



for run in RUNlist:
    print run
    filepkl = PKLdir + run + '/meanRSTbeforeagg.pkl'
    fid = open(filepkl,'r')
    LISTb = pickle.load(fid)
    fid.close()
    times_b[run] = LISTb[0]
    values_b[run] = LISTb[1][masktype]

    filepkl = PKLdir + run + '/meanRSTafteragg.pkl'
    fid = open(filepkl,'r')
    LISTa = pickle.load(fid)
    fid.close()
    times_a[run] = LISTa[0]
    values_a[run] = LISTa[1][masktype]

    for idate,dd in enumerate(times_b[run]):
#        times_b[run][idate] = dd.replace(hour=11,minute=59)
        times_b[run][idate] = dd-datetime.timedelta(minutes=1)

    tt[run] = times_a[run] + times_b[run]
    rr = np.array(tt[run])
    indt[run] = np.argsort(rr)

    tt[run].sort()
    

DICTalign = {
    'TRANSITION': 1,
    'RUN_TEST': -1,
    'RUN_2017Dard': -1,
    'everywhere': 1,
    'open': -1,
}


#import sys
#sys.exit(0)



for isub,sub in enumerate(OGS.P):
    print sub.name
    plt.close('all')
    fig,axs = plt.subplots(4,1,sharex=True,figsize=[12,10])
    

    plt.sca(axs[0])
    for run in RUNlist:
        valseries = np.concatenate((values_a[run][isub,:],values_b[run][isub,:]))
        val2plot = []
        for ii in range(len(tt[run])):
            val2plot.append(valseries[indt[run][ii]])
        plt.plot(tt[run],val2plot,label=run)
    
    plt.plot(sattimes,satvalues[isub,:],'.g',label='SAT')

    plt.ylabel('RST surface chl')
    plt.legend(loc='best')
    plt.title('Means over ' + sub.name + ' ' + masktype.upper())
    plt.grid()

    ax2 = axs[0].twinx()
    axs[0].set_zorder(ax2.get_zorder()+1)
    axs[0].patch.set_visible(False)

    plt.sca(ax2)
    plt.bar(sattimes,satcov[isub,:]/npoints[isub]*100, \
                width=2,align='center',color='lightgrey')
    plt.ylabel('Sat covergae %',color='grey')
    plt.ylim(0,100)
    # plt.grid(color='grey')


    plt.sca(axs[1])
    for run in RUNlist:
        plt.bar(times_b[run],values_a[run][isub,:]-values_b[run][isub,:], \
                width=2*DICTalign[run],align='edge',label=run)
    plt.ylabel('chl increments')
    plt.grid()
    # plt.legend(loc='best')

    plt.sca(axs[2])
    for run in RUNlist:
        plt.bar(misftimes[run],misfvalues[run][isub,:], \
                width=2*DICTalign[run],align='edge',label=run)
    plt.ylabel('innovations')
    plt.grid()
    # plt.legend(loc='best')

    plt.sca(axs[3])
    for run in RUNlist:
        plt.bar(misftimes[run],RMSvalues[run][isub,:], \
                width=2*DICTalign[run],align='edge',label=run)
    plt.ylabel('RMSbefore')
    plt.grid()
    plt.legend(loc='best')

    plt.tight_layout()

    plt.savefig(OUTdir + '/daTimes_2017Dard' + sub.name + masktype + '.png')





plt.show(block=False)




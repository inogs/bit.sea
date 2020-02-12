import numpy as np
from commons.submask import SubMask
from basins import V2 as OGS
import pickle
import datetime
import matplotlib.pyplot as plt


RUNlist = ['RUN_TEST','TRANSITION']
RUNlist = ['RUN_2017Dard','RUN_2017Dard_04']
masktype = 'open'
masktype = 'everywhere'

PKLdir = '/gpfs/scratch/userexternal/ateruzzi/ELAB_CFR_24/RST_DA/'
OUTdir = '/gpfs/scratch/userexternal/ateruzzi/ELAB_CFR_24/FIGURES_DA/'

Nsub = len(OGS.P.basin_list)


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



DICTalign = {
    'TRANSITION': 1,
    'RUN_TEST': -1,
    'RUN_2017Dard': -1,
    'RUN_2017Dard_04': 1,
    'everywhere': 1,
    'open': -1,
}


#import sys
#sys.exit(0)



for isub,sub in enumerate(OGS.P):
    print sub.name
    plt.close('all')
    fig,axs = plt.subplots(1,1,sharex=True,figsize=[12,7])
    
    for run in RUNlist:
        plt.bar(misftimes[run],RMSvalues[run][isub,:], \
                width=2*DICTalign[run],align='edge',label=run)
    plt.ylabel('RMSbefore')
    plt.grid()
    plt.legend(loc='best')

    plt.tight_layout()

    plt.savefig(OUTdir + '/misftimes_2017Dard' + sub.name + masktype + '.png')





plt.show(block=False)




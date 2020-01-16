import numpy as np
import matplotlib.pyplot as plt

from basins import V2 as OGS

INDIR = ''
OUTDIR = ''
INDIR = '/gpfs/scratch/userexternal/ateruzzi/ELAB_RACOAST2017/DAdelta/'
OUTDIR = '/gpfs/scratch/userexternal/ateruzzi/ELAB_RACOAST2017/DAdelta/'

INDELTA = INDIR + 'MEANdelta/'
INMISF = INDIR + 'MEANmisf/'
INCORR = INDIR + 'MEANcorr/'



deltaTime = np.load(INDELTA + '/meandelta_Time.npy')
deltaDate = np.load(INDELTA + '/meandelta_Date.npy')
misfTime = np.load(INMISF + '/meanmisf_Time.npy')
misfDate = np.load(INMISF + '/meanmisf_Date.npy')
corrTime = np.load(INCORR + '/meancorr_Time.npy')
corrDate = np.load(INCORR + '/meancorr_Date.npy')


masktypeLIST = deltaTime[0].keys()
masktypeLIST.sort()


for isub,sub in enumerate(OGS.P):
    print sub.name
    deltaforplot = {}
    misfforplot = {}
    corrforplot = {}
    for masktype in masktypeLIST:
        deltaforplot[masktype] = []
        misfforplot[masktype] = []
        corrforplot[masktype] = []

    for idate,dateDA in enumerate(deltaDate):
        for masktype in masktypeLIST:
            deltaforplot[masktype].append(deltaTime[idate][masktype][isub])

    for idate,dateDA in enumerate(misfDate):
        for masktype in masktypeLIST:
            misfforplot[masktype].append(misfTime[idate][masktype][isub])

    for idate,dateDA in enumerate(corrDate):
        for masktype in masktypeLIST:
            corrforplot[masktype].append(corrTime[idate][masktype][isub])

    plt.close('all')

    fig,axs = plt.subplots(len(masktypeLIST),1,figsize=[10,9],sharex=True)
    plt.suptitle(sub.name)

    for imask,masktype in enumerate(masktypeLIST):
        plt.sca(axs[imask])
        plt.title(masktype)
        plt.plot(misfDate,misfforplot[masktype],'x',label='Innovation')
        plt.plot(deltaDate,deltaforplot[masktype],'x',label='Increment from RST diff')
        plt.plot(corrDate,corrforplot[masktype],'x',label='Incr from corr file')

        plt.legend(loc='best')
        plt.grid()

    plt.savefig(OUTDIR + 'FIGURES/deltaDA' + sub.name + '.png')

plt.show(block=False)




# timedelta = []
# massdelta = []
# yearmassdelta = {}

# for itime,filein in enumerate(TLdelta.filelist):
#     massb = np.loadtxt(filein) # kmolC
#     massdelta_molC = massb[0]*10**3
#     massdelta_gC = massdelta_molC*12
#     massdelta_gC_m2 = massdelta_gC/areamed_m2
#     massdelta.append(massdelta_gC_m2)
#     timedelta.append(TLdelta.Timelist[itime])
#     yy = TLdelta.Timelist[itime].year
#     yearmassdelta[yy] += massdelta_gC_m2


# misfDate = np.load(INMISF + '/meanmisf_Date.npy')
# misfMean = np.load(INMISF + '/meanmisf_Time.npy')




# plt.close('all')

# fog,ax1 = plt.subplots(figsize=(12,8))

# ax2 = ax1.twinx()

# ax1.plot(misfDate,misfMean,'.-',color='g',label='misfit')
# ax1.set_ylabel('chl misf [mg chl/m^3]')
# plt.legend(loc='upper left')

# ax2.plot(timedelta,massdelta,'-',label='DA delta')
# ax2.set_ylabel('DA delta [gC/m^2]')
# plt.legend(loc='upper right')

# plt.title('P_' + vargroup)
# plt.grid()

# plt.show(block=False)

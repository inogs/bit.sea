import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Similar to plot_timeseries_STD.py,
    it works for four PFTs validation together
    It generates two pictures, one with four pfts in four separate axes,
    and one with four pfts plotted together.

    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = ''' input dir with 4 pkl files'''
                                )

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' Output dir of png files'''
                                )

    return parser.parse_args()

args = argument()
import pickle
import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from basins import V2 as OGS
from commons.utils import addsep

class filereader():
    def __init__(self, filename):
        
        fid = open(filename,'rb')
        LIST = pickle.load(fid)
        fid.close()
        TIMES,_,_,MODEL_MEAN,SAT___MEAN,_,_,MODEL__STD,SAT____STD,CORR, NUMB = LIST
        self.TIMES = TIMES
        self.MODEL_MEAN = MODEL_MEAN
        self.SAT___MEAN  = SAT___MEAN
        self.MODEL__STD  = MODEL__STD
        self.SAT____STD  = SAT____STD
        


INPUTDIR=addsep(args.inputdir)
OUTDIR=addsep(args.outdir)

lamSAT = [ 412, 443, 490, 510, 555, 670]

rrs_MODEL_mean=np.zeros(len(lamSAT))
rrs_MODEL__std=np.zeros(len(lamSAT))

rrs_SAT___mean=np.zeros(len(lamSAT))
rrs_SAT____std=np.zeros(len(lamSAT))


MONTH_list=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

VARLIST_rrs= []

for lam in lamSAT:
    VARLIST_rrs.append('RRS'+str(lam))

MATRIX_LIST_rrs=[filereader(INPUTDIR + var  +'_open_sea.pkl') for var in VARLIST_rrs]

for isub,sub in enumerate(OGS.P):

    plt.close('all')
    fig,axs = plt.subplots(4,3, gridspec_kw = {'wspace':0.2, 'hspace':0.4}, sharex=True, sharey=True, dpi=100,figsize=(8,9))


    for month in range(12):

        index_row=int(month / 3)
        index_col=int(month % 3)

        ax=axs[index_row, index_col]
        ax.title.set_text(MONTH_list[month])

        ax.set_xlim([400,700])
        if index_row == 3:
            ax.set_xlabel(r'$\lambda~[nm]$')
        if index_col == 0:
            ax.set_ylabel(r'$Rrs~[st^{-1}]$')

        for ivar,var in enumerate(VARLIST_rrs):
           rrs_MODEL_mean[ivar]=MATRIX_LIST_rrs[ivar].MODEL_MEAN[month,isub]
           rrs_MODEL__std[ivar]=MATRIX_LIST_rrs[ivar].MODEL__STD[month,isub]
           rrs_SAT___mean[ivar]=MATRIX_LIST_rrs[ivar].SAT___MEAN[month,isub]
           rrs_SAT____std[ivar]=MATRIX_LIST_rrs[ivar].SAT____STD[month,isub]

        ax.plot(lamSAT,rrs_MODEL_mean,'o', color='r')
        ax.plot(lamSAT,rrs_MODEL_mean-rrs_MODEL__std,':', color='r',alpha=0.5)
        ax.plot(lamSAT,rrs_MODEL_mean+rrs_MODEL__std,':', color='r',alpha=0.5)
        ax.fill_between(lamSAT,rrs_MODEL_mean-rrs_MODEL__std,rrs_MODEL_mean+rrs_MODEL__std,color='r',alpha=0.3)

        ax.plot(lamSAT,rrs_SAT___mean,'o', color='g')
        ax.plot(lamSAT,rrs_SAT___mean-rrs_SAT____std,':', color='g',alpha=0.5)
        ax.plot(lamSAT,rrs_SAT___mean+rrs_SAT____std,':', color='g',alpha=0.5)
        ax.fill_between(lamSAT,rrs_SAT___mean-rrs_SAT____std,rrs_SAT___mean+rrs_SAT____std,color='g',alpha=0.3)


    outfile="%srrs_spectra_%s.png" %(OUTDIR,sub.name)
    fig.savefig(outfile)

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
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from timeseries.plot import read_pickle_file
from basins import V2 as OGS
from commons.utils import addsep

class filereader():
    def __init__(self, filename):
       
        data, TL = read_pickle_file(filename)
        self.TIMES = TL.Timelist
        self.MODEL_MEAN = data[:,:,1,0,0]
        self.MODEL__STD = data[:,:,1,0,1]
        


INPUTDIR=addsep(args.inputdir)
OUTDIR=addsep(args.outdir)

lamMOD = [ 250, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625,
           650, 675, 700, 725, 775, 850, 950, 1050, 1150, 1250, 1350, 1450, 1550,
          1650, 1750, 1900, 2200, 2900, 3700 ];

Ed_mean=np.zeros(len(lamMOD))
Es_mean=np.zeros(len(lamMOD))

Ed__std=np.zeros(len(lamMOD))
Es__std=np.zeros(len(lamMOD))

MONTH_list=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

VARLIST_ed= []
VARLIST_es= []

for lam in lamMOD:
    VARLIST_ed.append('Ed_'+str(lam).zfill(4))
    VARLIST_es.append('Es_'+str(lam).zfill(4))

MATRIX_LIST_ed=[filereader(INPUTDIR + var  +'.pkl') for var in VARLIST_ed]
MATRIX_LIST_es=[filereader(INPUTDIR + var  +'.pkl') for var in VARLIST_es]


for isub,sub in enumerate(OGS.P):

    plt.close('all')
    fig,axs = plt.subplots(4,3, gridspec_kw = {'wspace':1.00, 'hspace':1.0})

    for month in range(12):

        index_row=int(month / 3)
        index_col=int(month % 3)

        ax=axs[index_row, index_col]
        ax.title.set_text(MONTH_list[month])

        ax.set_xlim([400,700])
        ax.set_xlabel(r'$\lambda$')
        ax.set_ylabel(r'$W~m^{-2}$')

        for ivar,var in enumerate(VARLIST_ed):
           Ed_mean[ivar]=MATRIX_LIST_ed[ivar].MODEL_MEAN[month,isub]
           Ed__std[ivar]=MATRIX_LIST_ed[ivar].MODEL__STD[month,isub]

        for ivar,var in enumerate(VARLIST_es):
           Es_mean[ivar]=MATRIX_LIST_es[ivar].MODEL_MEAN[month,isub]
           Es__std[ivar]=MATRIX_LIST_ed[ivar].MODEL__STD[month,isub]

        ax.plot(lamMOD,Ed_mean,'o', color='r')
        ax.plot(lamMOD,Ed_mean-Ed__std,':', color='r',alpha=0.5)
        ax.plot(lamMOD,Ed_mean+Ed__std,':', color='r',alpha=0.5)
        ax.fill_between(lamMOD,Ed_mean-Ed__std,Ed_mean+Ed__std,color='r',alpha=0.3)

        ax.plot(lamMOD,Es_mean,'o', color='b')
        ax.plot(lamMOD,Es_mean-Es__std,':', color='b',alpha=0.5)
        ax.plot(lamMOD,Es_mean+Es__std,':', color='b',alpha=0.5)
        ax.fill_between(lamMOD,Es_mean-Es__std,Es_mean+Es__std,color='b',alpha=0.3)


    outfile="%sed_es_spectra_%s.png" %(OUTDIR,sub.name)
    fig.savefig(outfile)

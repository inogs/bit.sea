# Generates images (in parallel) to compare different runs,
# using STAT_PROFILES directories

# Edit the USER SETTINS section below before launch.


import numpy as np
import pylab as pl
from basins import V2
from commons.mask import Mask
from timeseries.plot import read_pickle_file
import os

BFMv5_dict={     'Ac':'ALK',
                'ppn': 'netPPYc',
                'ppg':'ruPPYc',
                'ppb':'ruPBAc',
      'CaCO3flux_dic':'rcalCARc' }

try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False


class plot_container():
    def __init__(self, labelstring, color, path,maskobj):
        self.name=labelstring
        self.path= path
        self.mask = maskobj
        self.plotargs=color

    def load(self, varname):
        filename=self.path + varname + ".pkl"
        if varname in BFMv5_dict.keys():
            bfmv5_filename = self.path + BFMv5_dict[varname] + ".pkl"
            if (~os.path.exists(filename) & os.path.exists(bfmv5_filename)):
                filename=bfmv5_filename

        data, TL = read_pickle_file(filename)
        self.timelist=TL.Timelist
        self.values=data
        #              tim| sub basin|   coast| depth|  stat|
        

    def plot(self, axes_list,LEVELS,iSub):
        '''
        Arguments:
        * axes_list * a list of axis objects, the output of figure_generator.gen_structure()
        * LEVELS    * a list of floats, indicating in meters the depth to be plotted for each axes;
                      len(axes_list) must be equal to len(LEVELS) + 1
        * iSub      * integer, index of subbasin in STAT_PROFILES
        Returns : nothing
        '''
        iCoast=1 # open sea
        for iax, ax in enumerate(axes[:-1]):
            idepth = self.mask.getDepthIndex(LEVELS[iax])
            y = self.values[:,iSub,iCoast, idepth,0]
            if ~np.isnan(y).all():
                ax.plot(self.timelist,y,self.plotargs, label=self.name)
        ax = axes[-1] # the profile
        ax.plot(self.values[:,iSub,iCoast,:,0].mean(axis=0), self.mask.zlevels, self.plotargs, label=self.name )
        
class figure_generator():
    def __init__(self):
        self.LEFT_SIDE_AXES=[]
        self.ax_p=None

    def gen_structure(self, var, subname):
        '''
        Generates a figure structure
        Arguments :
        * var     * string
        * subname * string
        They are used just for the profile axes, at the left of the figure.


        Returns:
        * fig * figure object
        * AXES * a list of axes objects
        '''
        fig = pl.figure(figsize=(10, 10))
        ax1 = fig.add_subplot(421)
        ax2 = fig.add_subplot(423)
        ax3 = fig.add_subplot(425)
        ax4 = fig.add_subplot(427)
        
        self.LEFT_SIDE_AXES=[ax1, ax2, ax3, ax4]
        for iax, ax in enumerate(self.LEFT_SIDE_AXES):
            if iax<len(self.LEFT_SIDE_AXES)-1:
                ax.set_xticks([])
        
        ax_p = fig.add_subplot(122)
        title = "%s %s" %(var, subname)
        ax_p.set_title(title)
        ax_p.set_ylim([0, 1000])
        ax_p.invert_yaxis()
        self.ax_p = ax_p
        
        AXES=self.LEFT_SIDE_AXES[:]
        AXES.append(ax_p)
        return fig, AXES

    def add_text(self, LEVELS):
        '''
        Adds a text with info about depth
        In order to save space, the text is inside the figure instead of putting it on the title
        Arguments :
        * LEVELS    * a list of floats, indicating in meters the depth to be plotted for each axes;
                      len(axes_list) must be equal to len(LEVELS) + 1

        Returns : nothing
        '''
        for iax, ax in enumerate(self.LEFT_SIDE_AXES):
            Left,Right=ax.get_xlim()
            Bottom, Top =ax.get_ylim()
            title="depth = %dm" %LEVELS[iax]
            x=Left + (Right-Left)*0.03
            y=Top - (Top-Bottom)*0.05
            ax.text(x, y, title,fontsize=10,ha='left',va='top')
    def add_legend(self):
        self.ax_p.legend()
        
                    
##### USER SETTINGS #######################################
LOC="/marconi_work/OGS_dev_0/"


PATH6 = LOC + "DEGRADATION_4_70/CFR_PREVIOUS_RUNS/HC16/STAT_PROFILES/"
PATH8 = LOC + "DEGRADATION_4_70/TEST_08/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/"
PATH14= LOC + "DEGRADATION_4_70/TEST_14/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/"
PATH16= LOC + "DEGRADATION_4_70/TEST_16/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/"
PATH17= LOC + "DEGRADATION_4_70/TEST_17/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/"

Mask_4=Mask(LOC + "DEGRADATION_4_70/PREPROC/preproc_meshgen_forcings/mesh_gen/meshmask_470.nc",loadtmask=False)
Mask16=Mask(LOC + "DEGRADATION_4_70/POSTPROC/MASKS/meshmask16.nc",loadtmask=False)
Mask24=Mask(LOC + "DEGRADATION_4_70/POSTPROC/MASKS/meshmask24.nc",loadtmask=False)
OUTDIR= LOC + "DEGRADATION_4_70/POSTPROC/IMG/"

LEVELS=[0,50,100,150] #m

P14= plot_container('HC16_4_BFMv2', "r:"   ,PATH14, Mask_4)
P16= plot_container('HC16_4_bfmv5', "g"  , PATH16, Mask_4)
P17= plot_container('HC16_4_bfmv5_day_night', "k"  , PATH16, Mask_4)

PLOT_LIST=[P14,P16,P17]

VARLIST=["Ac","N1p", "N3n", "O2o", "P_l", "P_c", "DIC", "ppn", "P_n", "P_p", 'pH', "pCO2", "ppg", "ppb", "resPBAc", "CaCO3flux_dic"]
VARLIST_only_v5=['exPPYcR1', 'exPPYcR2','exPPYcR3','exPPYcR6',"resPPYc"]

##################################################################

for var in VARLIST[rank::nranks] :
    
    for p in PLOT_LIST: p.load(var)

    for iSub, sub in enumerate(V2.P):
        # if sub.name != "ion2" : continue
        outfile = "%sMultirun_Profiles.%s.%s.png" %(OUTDIR,var,sub.name)
        print outfile

        FigureGenerator=figure_generator()
        fig, axes= FigureGenerator.gen_structure(var, sub.name)
        
        for p in PLOT_LIST: p.plot(axes, LEVELS, iSub)
        FigureGenerator.add_text(LEVELS)
        FigureGenerator.add_legend()
        
        #fig.show()

        fig.savefig(outfile)
        pl.close(fig)


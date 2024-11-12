# module load mpi4py/1.3.1--intelmpi--2017--binary
# srun -N 1 -n 15 -A IscrB_MED21BIO_1 --time=30:00  --mem=50gb --partition=gll_usr_prod --pty bash

# Generates images (in parallel) to compare different runs,
# using STAT_PROFILES directories

# Edit the USER SETTINGS section below before launch.


import numpy as np
import matplotlib.pyplot as pl
from bitsea.basins import V2
from bitsea.commons.mesh import Mesh
from bitsea.timeseries.plot import read_pickle_file
import os

BFMv5_dict={     'Ac':'ALK',
                'ppn': 'netPPYc',
                'ppg':'ruPPYc',
                'ppb':'ruPBAc',
      'CaCO3flux_dic':'rcalCARc' }

BFMv2_dict={    'ALK':'Ac',
            'netPPYc': 'ppn' ,
           'netPPYc2': 'ppn' ,
             'ruPPYc': 'ppg',
             'ruPBAc':'ppb',
           'rcalCARc':'CaCO3flux_dic' }

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
        if varname in BFMv2_dict.keys():
            bfmv2_filename = self.path + BFMv2_dict[varname] + ".pkl"
            if (~os.path.exists(filename) & os.path.exists(bfmv2_filename)):
                filename=bfmv2_filename
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
                pl.gcf().autofmt_xdate()

        ax = axes[-1] # the profile
        ax.plot(self.values[:,iSub,iCoast,:,0].mean(axis=0), self.mask.zlevels, self.plotargs, label=self.name )
        
class figure_generator():
    def __init__(self):
        self.LEFT_SIDE_AXES=[]
        self.ax_p=None

    def gen_structure(self, var, subname,LEVELS):
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
                ax.set_xticklabels([])
            ax.grid()

        
        ax_p = fig.add_subplot(122)
        title = "%s %s" %(var, subname)
        ax_p.set_title(title)
        ax_p.set_ylim([0, 1000])
        y_ticklabels=[0,200,400,600,800,1000]
        y_ticklabels.extend(LEVELS)
        Y_TICK_LABELS=np.unique(y_ticklabels)
        ax_p.set_yticks(Y_TICK_LABELS)
        ax_p.grid()        
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
LOC="/g100_scratch/userexternal/gbolzon0/BI-HOURLY/"

PATH1 = LOC + "2H/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/"
PATH2 = LOC + "4H/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/"
PATH3 = LOC + "4H/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/"
PATH4 = LOC + "24H_orig/wrkdir/POSTPROC/output/AVE_FREQ_1/STAT_PROFILES/"

#Mask_4=Mesh.from_file(LOC + "DEGRADATION_4_70/PREPROC/preproc_meshgen_forcings/mesh_gen/meshmask_470.nc", read_e3t=True)
#Mask16=Mesh.from_file(LOC + "DEGRADATION_4_70/POSTPROC/MASKS/meshmask16.nc", read_e3t=True)
#Mask24=Mesh.from_file(LOC + "DEGRADATION_4_70/POSTPROC/MASKS/meshmask24.nc", read_e3t=True)

Mask24 = Mesh.from_file(LOC + "2H/wrkdir/MODEL/meshmask.nc", read_e3t=True)
OUTDIR= LOC + "multirun/spaghetti/offshore/"

LEVELS=[0,50,100,150] #m

P1= plot_container('2h', "r"   ,PATH1, Mask24)
P2= plot_container('4h', "b:"  ,PATH2, Mask24)
P3= plot_container('6h', "g"  , PATH3, Mask24)
P4= plot_container('24h', "k"  , PATH4, Mask24)


PLOT_LIST=[P1,P2,P3,P4]

VARLIST=["P_l","N3n","N1p","N4n","N5s","P_c","O2o"]
VARLIST=["N1p", "N3n", "N4n", "O2o", "P_l", "P_c", "DIC", "ppn", "ALK", "DIC", 'pH', "pCO2","O3c","O3h"]

# VARLIST=["DIC","pH","Ac","pCO2"]
# VARLIST_only_v5=['exPPYcR1', 'exPPYcR2','exPPYcR3','exPPYcR6',"resPPYc"]

##################################################################

for var in VARLIST[rank::nranks] :
    
    for p in PLOT_LIST: p.load(var)

    for iSub, sub in enumerate(V2.P):
        # if sub.name != "alb" : continue
        outfile = "%sMultirun_Profiles.%s.%s.png" %(OUTDIR,var,sub.name)
        print(outfile)

        FigureGenerator=figure_generator()
        fig, axes= FigureGenerator.gen_structure(var, sub.name,LEVELS)
        
        for p in PLOT_LIST: p.plot(axes, LEVELS, iSub)
        FigureGenerator.add_text(LEVELS)
        FigureGenerator.add_legend()
        
        # fig.show()
        # break

        fig.savefig(outfile)
        pl.close(fig)


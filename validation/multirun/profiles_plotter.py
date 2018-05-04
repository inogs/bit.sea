import pylab as pl
from basins import V2
from commons.mask import Mask
from timeseries.plot import read_pickle_file

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
    def __init__(self, labelstring, path,maskobj):
        self.name=labelstring
        self.path= path
        self.mask = maskobj

    def load(self, varname):
        filename=self.path + varname + ".pkl"
        data, TL = read_pickle_file(filename)
        self.timelist=TL.Timelist
        self.values=data
        #              tim| sub basin|   coast| depth|  stat|
        

    def plot(self, axes_list,LEVELS,iSub):        
        iCoast=1 # open sea
        for iax, ax in enumerate(axes[:-1]):
            idepth = self.mask.getDepthIndex(LEVELS[iax])
            ax.plot(self.timelist,self.values[:,iSub,iCoast, idepth,0],label=self.name)
        ax = axes[-1]
        ax.plot(self.values[:,iSub,iCoast,:,0].mean(axis=0), self.mask.zlevels, label=self.name )
        
class figure_generator():
    def __init__(self):
        pass
    

    def gen_structure(self, var, LEVELS, subname):
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

    def add_text(self):
        for iax, ax in enumerate(self.LEFT_SIDE_AXES):
            Left,Right=ax.get_xlim()
            Bottom, Top =ax.get_ylim()
            title="depth = %dm" %LEVELS[iax]
            x=Left + (Right-Left)*0.03
            y=Top - (Top-Bottom)*0.05
            ax.text(x, y, title,fontsize=10,ha='left',va='top')
    def add_legend(self):
        self.ax_p.legend()
        
                    


PATH1 = "/marconi_work/OGS_dev_0/DEGRADATION_4_70/TEST_04/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/"
PATH2 = "/marconi_work/OGS_dev_0/DEGRADATION_4_70/TEST_05/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/"
PATH3 = "/marconi_work/OGS_dev_0/DEGRADATION_4_70/TEST_06/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/"
PATH4 = '/marconi_scratch/userexternal/gbolzon0/TRANSITION_24/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/'
PATH5 = "/marconi_scratch/userexternal/gbolzon0/RA16_COAST/output/STAT_PROFILES/"
PATH6 = "/marconi_scratch/userexternal/gbolzon0/HC16/output/STAT_PROFILES/"

Mask_4=Mask("/marconi_work/OGS_dev_0/DEGRADATION_4_70/PREPROC/preproc_meshgen_forcings/mesh_gen/meshmask_470.nc",loadtmask=False)
Mask16=Mask("/marconi_work/OGS_dev_0/DEGRADATION_4_70/POSTPROC/MASKS/meshmask16.nc",loadtmask=False)
Mask24=Mask("/marconi_work/OGS_dev_0/DEGRADATION_4_70/POSTPROC/MASKS/meshmask24.nc",loadtmask=False)
OUTDIR="/marconi_work/OGS_dev_0/DEGRADATION_4_70/POSTPROC/IMG/"

LEVELS=[0,50,100,150] #m

P1 = plot_container('HC16_4_bfmv2'    , PATH1, Mask_4)
P2 = plot_container('HC16_4_bfmv5_st' , PATH2, Mask_4)
P3 = plot_container('HC16_4_bfmv5_ogs', PATH3, Mask_4)
P4 = plot_container('TRANS24'         , PATH4, Mask24)
P5 = plot_container('RA16'            , PATH5, Mask16)
P6 = plot_container('HC16'            , PATH6, Mask16)

PLOT_LIST=[P1,P2,P3,P4,P5,P6 ]

VARLIST=["N1p", "N3n", "O2o", "P_l", "P_c", "DIC"]

for var in VARLIST[rank::nranks] :
    
    for plot_obj in PLOT_LIST: plot_obj.load(var)

    for iSub, sub in enumerate(V2.P):
        FigureGenerator=figure_generator()
        fig, axes= FigureGenerator.gen_structure(var, LEVELS, sub.name)
        
        for plot_obj in PLOT_LIST: plot_obj.plot(axes, LEVELS, iSub)
        FigureGenerator.add_text()
        FigureGenerator.add_legend()
        

        #fig.show()

        outfile = "%sMultirun_Profiles.%s.%s.png" %(OUTDIR,var,sub.name)
        fig.savefig(outfile)
        print outfile
        pl.close(fig)     

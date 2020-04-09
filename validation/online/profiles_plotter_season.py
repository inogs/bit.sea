# module load mpi4py/1.3.1--intelmpi--2017--binary
# srun -N 1 -n 15 -A IscrB_MED21BIO_1 --time=30:00  --mem=50gb --partition=gll_usr_prod --pty bash

# Generates images (in parallel) to compare different runs,
# using STAT_PROFILES directories

# Edit the USER SETTINS section below before launch.
import argparse

def argument():
    parser = argparse.ArgumentParser(description='''
    plot post check and mis stats of BGC-Argo DA
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--indir', '-i',
                        type=str,
                        default=None,
                        required=True,
                        help="Dir with PKL dir")

    parser.add_argument('--outdir', '-o',
                        type=str,
                        default=None,
                        required=True,
                        help="FIGURES")

    parser.add_argument('--variables', '-v',
                        type=str,
                        default=None,
                        required=True,
                        help="FIGURES")

    parser.add_argument('--maskfile', '-m',
                        type=str,
                        default=None,
                        required=True,
                        help="meshmask.nc")

    parser.add_argument('--daterun', '-d',
                        type=str,
                        default=None,
                        required=True,
                        help="meshmask.nc")

    return parser.parse_args()

args = argument()


import numpy as np
import matplotlib.pyplot as pl
from basins import V2
from commons.mask import Mask
from commons import season
from commons.Timelist import TimeList
from commons import timerequestors
from commons.utils import addsep
from timeseries.plot import read_pickle_file
from datetime import datetime as date
from datetime import timedelta as timedelta
from xml.dom import minidom
import os

INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)
filevars = args.variables
MaskFile = args.maskfile
daterun8 = args.daterun


xmldoc = minidom.parse(filevars)

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
        

    def plot(self, axes_list,LEVELS,iSub, datetoday):
        '''
        Arguments:
        * axes_list * a list of axis objects, the output of figure_generator.gen_structure()
        * LEVELS    * a list of floats, indicating in meters the depth to be plotted for each axes;
                      len(axes_list) must be equal to len(LEVELS) + 1
        * iSub      * integer, index of subbasin in STAT_PROFILES
        Returns : nothing
        '''
        iCoast=1 # open sea
        currentyear = datetoday.year
        seasonObj = season.season()
        seasonind = seasonObj.findseason(datetoday)
        for iax, ax in enumerate(axes_list[:-1]):
            idepth = self.mask.getDepthIndex(LEVELS[iax])
            y = self.values[:,iSub,iCoast, idepth,0]
            datescurrentY = []
            for tt in self.timelist:
                ddY = tt.replace(year=currentyear)
                datescurrentY.append(ddY)
            if ~np.isnan(y).all():
                # ax.plot(self.timelist,y,self.plotargs, label=self.name)
                ax.plot(datescurrentY,y,self.plotargs, label=self.name)
                pl.gcf().autofmt_xdate()

        ax = axes_list[-1] # the profile
        seasonstr = seasonObj.SEASON_LIST_NAME[seasonind]
        seasreq = timerequestors.Clim_season(seasonind,seasonObj)
        sind,_ = TimeList(self.timelist).select(seasreq)
        maskseas = np.zeros_like(y,dtype=np.bool)
        maskseas[sind] = True
        ax.plot(self.values[maskseas,iSub,iCoast,:,0].mean(axis=0), \
                self.mask.zlevels, self.plotargs, \
                label=self.name + ' ' + seasonstr)
        
class figure_generator():
    def __init__(self):
        self.LEFT_SIDE_AXES=[]
        self.ax_p=None
        self.ax1=None

    def gen_structure(self, var, subname,LEVELS,datetoday):
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
        
        ax_p = fig.add_subplot(122)
        datestr = datetoday.strftime('%Y-%m-%d')
        title = "%s %s Last daterun %s" %(var,subname,datestr)
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
            ax.grid()

    def add_legend(self):
        self.ax_p.legend(loc='lower right')
        
    def add_seasonarea(self,datetoday):
        currenty = datetoday.year
        seasonObj = season.season()
        seasonind = seasonObj.findseason(datetoday)
        seasdates,_ = seasonObj.get_season_dates(seasonind)
        seas_start = seasdates.start_time.replace(year=currenty)
        seas_end   = seasdates.end_time.replace(year=currenty)

        datelastDA = datetoday-timedelta(days=2)

        for iax, ax in enumerate(self.LEFT_SIDE_AXES):
            ax.axvspan(seas_start,seas_end,facecolor='grey', alpha=0.2)
            ax.axvline(x=datelastDA,color='b',ls='--',alpha=0.5,lw=2)

                    
##### USER SETTINGS #######################################
LOC = INDIR

datetoday = date.now()
datetoday = date.strptime('20200415','%Y%m%d')
datetoday = date.strptime(daterun8,'%Y%m%d')

year = datetoday.year
LISTyear = np.arange(year-3,year+1)

PATH = {}
for yy in LISTyear:
    PATH[yy] =  LOC  + '/PKL' + np.str(yy) + '/'


Mask24=Mask(MaskFile)

LEVELS=[0,50,100,150] #m


PLOT_LIST = []
StyleLIST = [':g',':y',':c','-b']
for ii,yy in enumerate(LISTyear):
    PLOT_LIST.append(plot_container(np.str(yy),StyleLIST[ii],PATH[yy],Mask24))


VARLIST=["P_l"]

NN = xmldoc.getElementsByTagName("LayersMaps")
MM = NN[0].getElementsByTagName("plots")
VARLIST = []
for n in MM:
      VARLIST.append(str(n.attributes['var'].value))



##################################################################

for var in VARLIST[rank::nranks] :
#for var in [VARLIST[0]] :
    
    for p in PLOT_LIST: p.load(var)

    for iSub, sub in enumerate(V2.P):
        outfile = "%sMultirun_Profiles.%s.%s.png" %(OUTDIR,var,sub.name)
        print outfile

        FigureGenerator=figure_generator()
        fig, axes= FigureGenerator.gen_structure(var, sub.name,LEVELS,datetoday)
        
        #for p in PLOT_LIST: p.plot_timeser(axes, LEVELS, iSub, year)
        for p in PLOT_LIST: p.plot(axes, LEVELS, iSub, datetoday)
        FigureGenerator.add_text(LEVELS)
        FigureGenerator.add_legend()
        FigureGenerator.add_seasonarea(datetoday)
        
        # fig.show()
        # break

        fig.savefig(outfile)
        pl.close(fig)


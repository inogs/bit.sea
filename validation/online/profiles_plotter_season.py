import argparse

def argument():
    parser = argparse.ArgumentParser(description='''
    Generates spaghetti plots
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--inputdir', '-i',
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
                        help="Plotlist_bio.xml")

    parser.add_argument('--maskfile', '-m',
                        type=str,
                        default=None,
                        required=True,
                        help="meshmask.nc")

    parser.add_argument('--rundate', '-d',
                        type=str,
                        default=None,
                        required=True,
                        help="date in yyyymmdd format")

    return parser.parse_args()

args = argument()

import matplotlib
matplotlib.use("Agg")
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

INPUTDIR = addsep(args.inputdir)
OUTDIR = addsep(args.outdir)

xmldoc = minidom.parse(args.variables)


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
                if (tt.month,tt.day) != (2,29):  
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
        * LEVELS  * list of float or integers
        * datetoday * a datetime object
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
        title = "%s %s rundate %s" %(var,subname,datestr)
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
datetoday = date.strptime(args.rundate,'%Y%m%d')

year = datetoday.year
LISTyear = np.arange(year-2,year+1)

PATH = {}
for yy in LISTyear:
    PATH[yy] =  INPUTDIR  + '/PKL' + np.str(yy) + '/'


Mask24=Mask(args.maskfile)

LEVELS=[0,50,100,150] #m


PLOT_LIST = []
StyleLIST = [':g',':y',':c','-b']
for ii,yy in enumerate(LISTyear):
    PLOT_LIST.append(plot_container(np.str(yy),StyleLIST[ii],PATH[yy],Mask24))


NN = xmldoc.getElementsByTagName("LayersMaps")
MM = NN[0].getElementsByTagName("plots")
VARLIST = []
for n in MM:
    VARLIST.append(str(n.attributes['var'].value))


##################################################################

for var in VARLIST[rank::nranks] :
    
    for p in PLOT_LIST: p.load(var)

    for iSub, sub in enumerate(V2.P):
        outfile = "%sProfiles.%s.%s.png" %(OUTDIR,var,sub.name)
        print outfile

        FigureGenerator=figure_generator()
        fig, axes= FigureGenerator.gen_structure(var, sub.name,LEVELS,datetoday)
        
        for p in PLOT_LIST: p.plot(axes, LEVELS, iSub, datetoday)
        FigureGenerator.add_text(LEVELS)
        FigureGenerator.add_legend()
        FigureGenerator.add_seasonarea(datetoday)


        fig.savefig(outfile)
        pl.close(fig)


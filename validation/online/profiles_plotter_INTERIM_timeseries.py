import argparse

import argparse

def argument():
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = 'The directory that contains pkl files')    

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = 'The directory which will contain spaghetti plot files')

    parser.add_argument(   '--timelabel', '-t',
                                type = str,
                                required =True,
                                help = '''
                                name of last month that will be plotted alone
                                EX: December 2022
                                ''')
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
import pylab as pl
import matplotlib.cm as cm
import matplotlib.dates as mdates

INPUTDIR = addsep(args.inputdir)
OUTPUTDIR = addsep(args.outdir)
time_label = addsep(args.timelabel)

VARLIST=["ALK","DIC","exR2ac","N1p","N3n","N4n","N5s","O2o","P_c","pCO2","pH","P_l","ppn",'Z_c']
#xmldoc = minidom.parse(args.variables)


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
    def __init__(self, labelstring, color, path,maskobj,label):
        self.name=labelstring
        self.path= path
        self.mask = maskobj
        self.plotargs=color
        self.label = label

    def load(self, varname):
        filename=self.path + varname + ".pkl"
        data, TL = read_pickle_file(filename)
        self.timelist=TL.Timelist
        self.values=data
        #              tim| sub basin|   coast| depth|  stat|
        

    def plot(self, axes_list,LEVELS,iSub,datetoday):
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
                ax.plot(datescurrentY,y,color=self.plotargs)#, label=self.name)
                pl.gcf().autofmt_xdate()

        ax = axes_list[-1] # the profile
        seasonstr = seasonObj.SEASON_LIST_NAME[seasonind]
        seasreq = timerequestors.Clim_season(seasonind,seasonObj)
        sind,_ = TimeList(self.timelist).select(seasreq)
        maskseas = np.ones_like(y,dtype=np.bool)
        #maskseas[sind] = True
        if self.label is not None:
            ax.plot(self.values[maskseas,iSub,iCoast,:,0].mean(axis=0),self.mask.zlevels, color=self.plotargs,label=self.label)
            #ax.plot(self.values[-1,iSub,iCoast,:,0],self.mask.zlevels, color='g',ls='--',label=self.label)
        else:
            ax.plot(self.values[maskseas,iSub,iCoast,:,0].mean(axis=0),self.mask.zlevels, color=self.plotargs)
       
    def plot_extra(self,axes_list,LEVELS,iSub):
        iCoast=1 #Open sea

        axes_list[-1].plot(self.values[-1,iSub,iCoast,:,0],self.mask.zlevels, color='b',ls='--',label=time_label)

class figure_generator():
    def __init__(self):
        self.LEFT_SIDE_AXES=[]
        self.ax_p=None
        self.ax1=None

    def gen_structure(self, var, subname,LEVELS,datetoday):
        '''
        Generates a figure structur)#, label=self.name)e
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
        title = "%s %s" %(var,subname)#,datestr)
        ax_p.set_title(title)
        ax_p.set_ylim([0, 1000])
        y_ticklabels=[0,200,400,600,800,1000]
        y_ticklabels.extend(LEVELS)  #aggiunge livelli se sono passati
        Y_TICK_LABELS=np.unique(y_ticklabels)
        ax_p.set_yticks(Y_TICK_LABELS)
        ax_p.grid()        
        ax_p.invert_yaxis()
        self.ax_p = ax_p
        
        AXES=self.LEFT_SIDE_AXES[:]
        AXES.append(ax_p)
        ax4.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
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

    #def add_season(self):

     #       self.ax_p.plot(self.values[-1,iSub,iCoast,:,0],self.mask.zlevels, color='g',ls='--',label=self.label)
                    
##### USER SETTINGS #######################################
datetoday = date.strptime('20211001','%Y%m%d')

#year = datetoday.year
#LISTyear = np.arange(year-2,year+1)

LISTyear= range(1999,2018)

PATH = {}
for yy in LISTyear :
    PATH[yy] =  INPUTDIR  + '/' + np.str(yy) +'/'
PATH_2019=INPUTDIR + '/' + np.str(2019) +'/'
PATH_2020=INPUTDIR + '/' + np.str(2020) +'/'
PATH_2021=INPUTDIR + '/' + np.str(2021) +'/'
PATH_2022=INPUTDIR + '/' + np.str(2022) +'/'

Mask24=Mask('/g100_scratch/userexternal/gcoidess/REANALISI_INTERIM/wrkdir/MODEL/meshmask.nc')

LEVELS=[0,50,100,150] #m
ncolors=19

COLOR=cm.gray_r(np.linspace(0.2,0.3,ncolors))

PLOT_LIST = []

for ii,yy in enumerate(LISTyear):
    PLOT_LIST.append(plot_container(np.str(yy),COLOR[ii],PATH[yy],Mask24,None))

PLOT_LIST.append(plot_container('2019',COLOR[18],PATH_2019,Mask24,'1999-2019'))
PLOT_LIST.append(plot_container('2020','r',PATH_2020,Mask24,'2020'))
PLOT_LIST.append(plot_container('2021','g',PATH_2021,Mask24,'2021'))
PLOT_LIST.append(plot_container('2022','b',PATH_2022,Mask24,'2022'))

#NN = xmldoc.getElementsByTagName("LayersMaps")
#MM = NN[0].getElementsByTagName("plots")
#VARLIST = []
#for n in MM:
#    VARLIST.append(str(n.attributes['var'].value))


##################################################################

for var in VARLIST[rank::nranks] :
    
    for p in PLOT_LIST: p.load(var)
    for iSub, sub in enumerate(V2.P):
            outfile = "%sProfiles.%s.%s.png" %(OUTDIR,var,sub.name)
            print (outfile)

            FigureGenerator=figure_generator()
            fig, axes= FigureGenerator.gen_structure(var, sub.name,LEVELS,datetoday)
            import sys
        
            for ip,p in enumerate(PLOT_LIST):
                print(ip)
                p.plot(axes, LEVELS, iSub,datetoday)
            PLOT_LIST[-1].plot_extra(axes,LEVELS,iSub)
        
            FigureGenerator.add_text(LEVELS)
            FigureGenerator.add_legend()
            #sys.exit()
            #FigureGenerator.add_seasonarea(datetoday)
            #FigureGenerator.add_season()


            fig.savefig(outfile)
            pl.close(fig)
            #sys.exit()

import argparse
from datetime import date

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates float validation data on a single week, from one single chain run.
    Week is centered on tuesday.
    Produces a single file, containing bias, rmse and number of measurements for each subbasin and layer
    for chlorophyll, nitrate and oxygen.
    In this approach we define the measurement as mean on layer of the float profile values.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--pathV4','-V4',
                                type = str,
                                required = True,
                                help = 'Directory containing STAT_PROFILES/')

    parser.add_argument(   '--pathV5','-V5',
                                type = str,
                                required = True,
                                help = 'Directory containing STAT_PROFILES/')

    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc",
                                required = False,
                                help = ''' Path of maskfile''')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

import matplotlib
#matplotlib.use('Agg')
import matplotlib.dates as mdates
import numpy as np
import pylab as pl
from basins import V2
from commons.mask import Mask
from timeseries.plot import read_pickle_file
from dateutil.relativedelta import relativedelta
from datetime import datetime
import os
from commons.time_interval import TimeInterval
from commons import timerequestors
from commons.utils import addsep
from commons import season
Season_obj=season.season()

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
        self.timelist=TL
        self.values=data
        #              tim| sub basin|   coast| depth|  stat|
        

    def plot_timeseries(self, axes_list,LEVELS,iSub):
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
                ax.plot(self.timelist.Timelist,y,self.plotargs, label=self.name)

        
    def plot_mean_profiles(self, axes_list,LEVELS,iSub):
        iCoast=1 # open sea        
        ax = axes[-1] # the profile
        color=(0.7,0.7, 0.7)
        nFrames=self.timelist.nTimes
        for k in range(7):
            ax.plot(self.values[nFrames-k-1,iSub,iCoast,:,0], self.mask.zlevels, color=color, linestyle='--')
        ax.plot(self.values[nFrames-k-1,iSub,iCoast,:,0], self.mask.zlevels, color=color, linestyle='--', label="last 7 days")
        
        
        seas_colors=['k','g','m','c']
        for iSeas in range(4):
            seas_req=timerequestors.Clim_season(iSeas, Season_obj)
            ii_seas,w = self.timelist.select(seas_req)
            ax.plot(self.values[ii_seas,iSub,iCoast,:,0].mean(axis=0), self.mask.zlevels, seas_colors[iSeas], label="clim " + seas_req.string)

        
        
class figure_generator():
    def __init__(self):
        self.LEFT_SIDE_AXES=[]
        self.ax_p=None

    def gen_structure(self, TI, var, subname,LEVELS):
        '''
        Generates a figure structure
        Arguments :
        * TI      * a TimeInterval Object
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
            ax.set_xlim(TI.start_time, TI.end_time)
            if iax<len(self.LEFT_SIDE_AXES)-1:
                ax.set_xticks([])
        ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,4,7,10]))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
        xlabels = ax.get_xticklabels()
        pl.setp(xlabels, rotation=30)
        
        ax_p = fig.add_subplot(122)
        title = "%s %s" %(var, subname)
        ax_p.set_title(title)
        ax_p.set_ylim([0, 600])
        y_ticklabels=[0,200,400,600]
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
        self.LEFT_SIDE_AXES[0].legend()
        
                    


PATHV4 = addsep(args.pathV4)
PATHV5 = addsep(args.pathV5)
TheMask=Mask(args.maskfile,loadtmask=False)
OUTDIR= addsep(args.outdir)

LEVELS=[0,50,100,150] #m

P4= plot_container('V4', "b"   ,PATHV4, TheMask)
P5= plot_container('V5', "r"   ,PATHV5, TheMask)
PLOT_LIST=[P4, P5]

VARLIST=["N1p", "N3n", "O2o", "P_l", "P_c", "N4n", "N5s", "DIC", "Ac",  "ppn",  'pH', "pCO2", "CO2airflux"]

dateend=datetime.strptime("20190415",'%Y%m%d')
Graphic_DeltaT = relativedelta(months=18)
datestart = dateend -Graphic_DeltaT
TI =TimeInterval.fromdatetimes(datestart, dateend)
##################################################################

for var in VARLIST[rank::nranks] :
    
    for p in PLOT_LIST: p.load(var)

    for iSub, sub in enumerate(V2.P):
        # if sub.name != "ion2" : continue
        outfile = "%sMultirun_Profiles.%s.%s.png" %(OUTDIR,var,sub.name)
        print outfile

        FigureGenerator=figure_generator()
        fig, axes= FigureGenerator.gen_structure(TI, var, sub.name,LEVELS)
        
        P4.plot_timeseries(axes,LEVELS,iSub)
        P5.plot_timeseries(axes,LEVELS,iSub)
        P5.plot_mean_profiles(axes, LEVELS, iSub)
        #for p in PLOT_LIST: p.plot(axes, LEVELS, iSub)
        
        
        FigureGenerator.add_text(LEVELS)
        FigureGenerator.add_legend()

        fig.savefig(outfile)
        pl.close(fig)



import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates spaghetti plots, by using a 

    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc",
                                required = True,
                                help = ''' Path of maskfile''')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")
    parser.add_argument(   '--settings_file', '-f',
                                type = str,
                                default = "profiles_plotter_user_settings.txt",
                                required = False,
                                help = "")    

    return parser.parse_args()

args = argument()

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as pl
from bitsea.basins import V2
from bitsea.commons.mesh import Mesh
from bitsea.timeseries.plot import read_pickle_file
import os
from bitsea.commons.utils import addsep

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

Mask24 = Mesh(args.maskfile, read_e3t=True)
OUTDIR= addsep(args.outdir)

mydtype=np.dtype([ ('label','U20'), ('linestyle','U5'), ('PATH','U1024')  ])

A=np.loadtxt(args.settings_file,dtype=mydtype)

PLOT_LIST = []
for p in A:
    P = plot_container(p['label'],p['linestyle'],addsep(p['PATH']),Mask24)
    PLOT_LIST.append(P)



LEVELS=[0,50,100,150] #m
VARLIST=["N1p", "N3n", "N4n", "O2o", "P_l", "P_c", "DIC", "ppn", "ALK", "DIC", 'pH', "pCO2","O3c","O3h",'CO2airflux']


##################################################################

for var in VARLIST:
    os.system('mkdir -p ' + OUTDIR + var +"/")
    for p in PLOT_LIST: p.load(var)

    for iSub, sub in enumerate(V2.P):
        # if sub.name != "alb" : continue
        outfile = "%s%s/Multirun_Profiles.%s.%s.png" %(OUTDIR,var,var,sub.name)
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


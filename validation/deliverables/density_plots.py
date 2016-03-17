import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates density plots
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                default = '',
                                help = 'Output Images directory')


    return parser.parse_args()


args = argument()

import pylab as pl
import scipy.io.netcdf as NC
import numpy as np
import os
import sys
from commons.time_interval import TimeInterval

from profiler_RA import *

import basins.OGS as OGS
from instruments.instruments import static_Selector
from commons.layer import Layer
M = Matchup_Manager(T_INT,INPUTDIR,BASEDIR)

maskfile    = os.getenv("MASKFILE");
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()

OUTPUTDIR=args.outdir
os.system('mkdir -p ' + OUTPUTDIR)

#modelvarname='N1p'
#modelvarname='N3n'
modelvarname='O2o'
for isub,sub in enumerate(OGS.P):
    print sub.name
    Profilelist=static_Selector(modelvarname,T_INT,sub)
    Matchup_basin = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=True)
# arg in () corrisponde al n. di bin dell'istogramma
    #fig,ax =  Matchup_basin.densityplot(10)
    fig,ax =  Matchup_basin.densityplot2(modelname='RAN',refname='REF',units='mmol P/m$^3$',sub=sub.name.upper())
    #fig,ax =  Matchup_basin.densityplot2(modelname='RAN',refname='REF',units='mmol N/m$^3$',sub=sub.name.upper())
    #fig,ax =  Matchup_basin.densityplot2(modelname='RAN',refname='REF',units='mmol O$_2$/m$^3$',sub=sub.name.upper())
    maxval=np.max(Matchup_basin.Ref.max(),Matchup_basin.Model.max())
    ax.set_xlim([0,maxval])
    ax.set_ylim([0,maxval])
    ax.set_xlabel('RAN data [mmol P/m$^3$]').set_fontsize(14)
    ax.set_ylabel('REF data [mmol P/m$^3$]').set_fontsize(14)
    ax.set_title(sub.name.upper() + ' - TOT n. matchups= ' + str(Matchup_basin.number()))
    ax.grid(True)
    fig.show()
    fig.savefig(OUTPUTDIR+'densplot_'+modelvarname+'_'+sub.name+'.png')
    sys.exit()

    #fig=pl.figure(num=None,dpi=100,facecolor='w',edgecolor='k')
    #ax=fig.add_subplot(1,1,1)
    #pl.plot(Matchup_basin.Ref,Matchup_basin.Depth,'.r',label='REF')
    #pl.plot(Matchup_basin.Model,Matchup_basin.Depth,'.k',label='RAN')
    #pl.gca().invert_yaxis()
    #ax.legend(loc="best",labelspacing=0, handletextpad=0,borderpad=0.1)
    #leg = pl.gca().get_legend()
    #ltext  = leg.get_texts()
    #ax.set_xlabel('[mmol P/m$^3$]').set_fontsize(14)
    #ax.set_xlabel('[mmol N/m$^3$]').set_fontsize(14)
    #ax.set_xlabel('[mmol O$_2$/m$^3$]').set_fontsize(14)
    #ax.set_ylabel('[m]').set_fontsize(14)
    #ax.set_title(sub.name.upper() + ' ' + modelvarname)
    #pl.show()
    #fig.savefig(OUTPUTDIR+'vertprof_'+modelvarname+'_'+sub.name+'.png')    

    print sub.name,'correlation= ',Matchup_basin.correlation()

    fig, ax = Matchup_basin.taylorplot(dpi=72)

    m2 = Matchup_basin.subset(Layer(0,50))



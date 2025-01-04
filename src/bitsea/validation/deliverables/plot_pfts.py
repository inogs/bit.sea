import argparse
from bitsea.utilities.argparse_types import date_from_str
from bitsea.utilities.argparse_types import existing_dir_path
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
                                type = existing_dir_path,
                                required = True,
                                help = ''' input dir with '''
                                )

    parser.add_argument(   '--outdir', '-o',
                                type = existing_dir_path,
                                required = True,
                                help = ''' Output dir of png files'''
                                )
    parser.add_argument(   '--datestart', '-s',
                                type = date_from_str,
                                required =True,
                                help = ''' Date start for time interval to consider for validation, format %Y%m%d'''
                                )
    parser.add_argument(   '--dateend', '-e',
                                type = date_from_str,
                                required =True,
                                help = ''' Date end for time interval to consider for validation,format %Y%m%d'''
                                )
    return parser.parse_args()

args = argument()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from bitsea.basins import V2 as OGS
from bitsea.validation.deliverables import netcdf_validation_file
from bitsea.commons.Timelist import TimeInterval

TI = TimeInterval.fromdatetimes(args.datestart, args.dateend)

INPUTDIR=args.inputdir
OUTDIR=args.outdir

VARLIST= ["P1l","P2l","P3l","P4l"]
PFT_NAME=['Diatoms','Nanophytoplankton','Picophytoplankton','Dinoflagellates']
COLOR  = ['tab:blue','tab:orange','tab:green','tab:purple']
LIGHTCOLOR  = ['lightsteelblue','moccasin','palegreen','plum']



MATRIX_LIST=[netcdf_validation_file.dir_reader(TI, INPUTDIR,var,"open_sea") for var in VARLIST]

#P1l = dir_reader(INPUTDIR,"P1l","open_sea")
#P2l = dir_reader(INPUTDIR,"P2l","open_sea")



for isub,sub in enumerate(OGS.P):
    pl.close('all')
    fig = pl.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(411)
    ax2 = fig.add_subplot(412)
    ax3 = fig.add_subplot(413)
    ax4 = fig.add_subplot(414)

    AXES_LIST=[ax1, ax2, ax3, ax4]


    model_label='MODEL'
    units="[mg/m$^3$]"
    var_label = 'Sat' + " " + units


    for ivar,var in enumerate(VARLIST):
        outfile = OUTDIR / f"pfts_separated_{sub.name}.png"
        color = COLOR[ivar]
        ax = AXES_LIST[ivar]
        Pl = MATRIX_LIST[ivar]
        ax.plot(Pl.TIMES,Pl.SAT___MEAN[:,isub],'o', color=color)
        ax.fill_between(Pl.TIMES,Pl.SAT___MEAN[:,isub]-Pl.SAT____STD[:,isub],Pl.SAT___MEAN[:,isub]+Pl.SAT____STD[:,isub],color=LIGHTCOLOR[ivar])
        ax.plot(Pl.TIMES,Pl.MODEL_MEAN[:,isub],'-',color=color)
        ax.plot(Pl.TIMES,Pl.MODEL_MEAN[:,isub]-Pl.MODEL__STD[:,isub],':',color=color)
        ax.plot(Pl.TIMES,Pl.MODEL_MEAN[:,isub]+Pl.MODEL__STD[:,isub],':',color=color)
        ax.grid(True)
        ax.set_ylabel("%s" %(PFT_NAME[ivar])).set_fontsize(14)
        ax.tick_params(axis='both', which='major', labelsize=13)
        ax.set_ylim(0,0.25)
        if ax != ax4: ax.set_xticklabels([])

    ax1.set_title("PFTs  %s - %s" %(sub.name.upper(), units ) ).set_fontsize(15)
    xlabels = ax4.get_xticklabels()
    #pl.setp(xlabels, rotation=30)

    fig.savefig(outfile)

    # plot with all pfts together
    outfile = OUTDIR / f"pfts_{sub.name}.png"

    fig = pl.figure(figsize=(10, 4))
    ax = fig.add_subplot(111)
    for ivar,var in enumerate(VARLIST):
        Pl = MATRIX_LIST[ivar]
        ax.plot(Pl.TIMES,Pl.MODEL_MEAN[:,isub],'-',color=COLOR[ivar], label=PFT_NAME[ivar])
        ax.plot(Pl.TIMES,Pl.SAT___MEAN[:,isub],'--', color=COLOR[ivar])


    ax.plot(Pl.TIMES,Pl.SAT___MEAN[:,isub],'--',color=COLOR[ivar])
    ax.set_ylabel("PFTs %s" %(units)).set_fontsize(15)
    ax.set_ylim(0,0.25)
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.legend()
    ax.grid(True )
    ax.set_title("PFTs %s" %(sub.name.upper()  ) ).set_fontsize(15)
    fig.savefig(outfile)


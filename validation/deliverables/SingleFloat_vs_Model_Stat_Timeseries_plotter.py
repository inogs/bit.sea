import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Needs a profiler.py, already executed.

    Produces 3 png files, containing timeseries for some statistics, for each wmo.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc",
                                required = False,
                                help = ''' Path of maskfile''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

import numpy as np
from commons.mask import Mask
from instruments import lovbio_float as bio_float
from instruments.matchup_manager import Matchup_Manager
import matplotlib.pyplot as plt
from commons.utils import addsep
from profiler import ALL_PROFILES,TL,BASEDIR
from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader

def fig_setup(wmo,Lon,Lat,var):
    from layer_integral import coastline

    fig = plt.figure()
    ax0 = plt.subplot2grid((4, 3), (0, 0), colspan=2)
    ax1 = plt.subplot2grid((4, 3), (0, 2))
    ax2 = plt.subplot2grid((4, 3), (1, 0), colspan=3)
    ax3 = plt.subplot2grid((4, 3), (2, 0), colspan=3)
    ax4 = plt.subplot2grid((4, 3), (3, 0), colspan=3)
    axs = [ax0, ax1, ax2, ax3, ax4]

    fig.set_size_inches(10,15)
    c_lon,c_lat=coastline.get()

#    list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo_list[j])
    ax0.plot(c_lon,c_lat,'k')
    ax0.plot(Lon,Lat,'r.')
    ax0.set_title("TRAJECTORY of FLOAT " + wmo , color = 'r')
#    ind_max_sup=plotmat[0,:].argmax()
    
#    print Lon[ind_max_sup],Lat[ind_max_sup]
#    ax0.plot(Lon[ind_max_sup],Lat[ind_max_sup],'g.')
#    ax0.plot(Lon[0],Lat[0],'bx')
    ax0.set_xlim([-10,36])
    ax0.set_ylabel("LAT",color = 'k')
    ax0.set_xlabel("LON",color = 'k')

    return fig, axs


TheMask=Mask(args.masfile)
INDIR = addsep(args.inputdir)
OUTDIR = addsep(args.outdir)

VARLIST = ['P_l','N3n','O2o']
nVar = len(VARLIST)
METRICS = ['Int_0-200','Corr','DCM','z_01','Nit_1']
nStat = len(METRICS)

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
wmo_list=bio_float.get_wmo_list(ALL_PROFILES)

izmax = TheMask.getDepthIndex(200) # Max Index for depth 200m

for wmo in wmo_list:
    INPUT_FILE = INDIR + wmo + ".nc"
    print INPUT_FILE
    A = ncreader(INPUT_FILE)
    wmo_track_list = bio_float.filter_by_wmo(ALL_PROFILES,wmo)
    nP = len(wmo_track_list)
    Lon = np.zeros((nP,), np.float64)
    Lat = np.zeros((nP,), np.float64)
    for ip, p in enumerate(wmo_track_list):
        Lon[ip] = p.lon
        Lat[ip] = p.lat
    times = [p.time for p in wmo_track_list]
    for var in VARLIST:
        OUTFILE = OUTDIR + var + "_" + wmo + ".png"
        print OUTFILE
        fig, axes = fig_setup(wmo,Lon,Lat,var)
        if (var == "P_l"): 
            model, ref =A.plotdata(var,'Int_0-200')
            axes[2].plot(times,  ref,'r')
            axes[2].plot(times,model,'b')
        fig.savefig(OUTFILE)

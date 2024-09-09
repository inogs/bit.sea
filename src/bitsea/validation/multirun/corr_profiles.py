import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Produces Hovmoeller png files, containing the float trajectory and the two CHLA Hovmoeller for float and model for each wmo.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/marconi_scratch/userexternal/lfeudale/Maskfiles/meshmask.nc",
                                required = False,
                                help = ''' Path of maskfile''')

    parser.add_argument(   '--maskfilephys', '-p',
                                type = str,
                                default = "/marconi_scratch/userexternal/lfeudale/Maskfiles/meshmask.nc",
                                required = False,
                                help = ''' Path of physics maskfile''')


    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    parser.add_argument(   '--depth', '-d',
                                type = float,
                                required = False,
                                default = 200,
                                help  = 'The level chosen to filter data' )

    return parser.parse_args()

args = argument()

import numpy as np
import scipy.io.netcdf as NC

from commons.layer import Layer
from commons.mask import Mask
from commons.utils import addsep
from metrics import *

from instruments import superfloat as bio_float
from instruments.var_conversions import FLOATVARS
from instruments.matchup_manager import Matchup_Manager

from profiler import *

from layer_integral import coastline
from basins import OGS

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from timeseries.plot import *
import numpy.ma as ma
from mhelpers.pgmean import PLGaussianMean
import seawater
import pickle

from SingleFloat_vs_Model_Stat_Timeseries_IOnc import ncreader

TheMask = Mask(args.maskfile)
TheMask_phys = Mask(args.maskfilephys)

OUTDIR = addsep(args.outdir)
depth_lim = args.depth

run = 'PHY_A'
#run = "PHY_S"


font_s =  15 
label_s = 15

def my_Hovmoeller_diagram(plotmat, xs,ys, fig=None, ax=None):
    if (fig is None) or (ax is None):
        fig , ax = plt.subplots()
    quadmesh = ax.pcolormesh(xs, ys, plotmat,shading='gouraud')# default is 'flat'
    #Inform matplotlib that the x axis is made by dates
    ax.xaxis_date()
    ax.invert_yaxis()
    return fig, ax, quadmesh

def readModelProfile(filename,var, wmo):
    ncIN = NC.netcdf_file(filename,'r')
    M = ncIN.variables[var].data.copy()
    iProfile = ncIN.CruiseIndex.rsplit(", ").index(wmo)
    ncIN.close()
    Profile = M[iProfile,:]
    return Profile

def get_level300(TheMask):

    i0 = TheMask.getDepthIndex(300)
    i1 = i0 + 1
      
    diff0 = 300 - TheMask.zlevels[i0]
    diff1 = 300 - TheMask.zlevels[i1]
    data = [(i0,diff0),(i1,diff1)]
    ix,datamin = min(data, key=lambda t: t[1])
    return ix

def get_level200(TheMask):

    i0 = TheMask.getDepthIndex(200)
    i1 = i0 + 1
      
    diff0 = 200 - TheMask.zlevels[i0]
    diff1 = 200 - TheMask.zlevels[i1]
    data = [(i0,diff0),(i1,diff1)]
    ix,datamin = min(data, key=lambda t: t[1])
    return ix

def get_level_depth(TheMask,dd=200):

    i0 = TheMask.getDepthIndex(dd)
    i1 = i0 + 1

    diff0 = dd - TheMask.zlevels[i0]
    diff1 = dd - TheMask.zlevels[i1]
    data = [(i0,diff0),(i1,diff1)]
    ix,datamin = min(data, key=lambda t: t[1])
    return ix


# max_depth = get_level300(TheMask)
# max_depthp = get_level300(TheMask_phys)
max_depth = get_level_depth(TheMask,dd=depth_lim)
max_depthp = get_level_depth(TheMask_phys,dd=depth_lim)

VARLIST_PHYS = ['votemper','vosaline']
VARLIST_BGC = ['N3n','Chla']
VARLIST_BGC = ['N3n','P_l']
Adj = {
	'P_l':  True,
	'Chla': True,
	'O2o':  False,
	'N3n':  True,
	'votemper':  False,
	'vosaline':  False,
}

dofloat = True
meanObj11 = PLGaussianMean(11,1.0)

T_start = DATESTART
T_end   = DATE__END
TI1 = T_INT

Profilelist_1=bio_float.FloatSelector(None,TI1,OGS.med)
wmo_list=bio_float.get_wmo_list(Profilelist_1)

MM = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

print wmo_list

for j,wmo in enumerate(wmo_list):
    print(wmo)

    list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo)
    nP = len(list_float_track)
    Lon = np.zeros((nP,), np.float64)
    Lat = np.zeros((nP,), np.float64)
    nlev= depth_lim/5 + 1
    NewPres_5m=np.linspace(0,depth_lim,nlev)
#    NewPres_5m=np.linspace(0,300,61)
#    NewPres_5m=np.linspace(0,200,41)
    nPnewpres = len(NewPres_5m)

    timelabel_list = list()

    corrfloat = {}
    corrmod = {}
    covnfloat = {}
    covnmod = {}
    RMSD = {}
    for varb in VARLIST_BGC: 
        RMSD[varb] = list()
    for varp in VARLIST_PHYS: 
        RMSD[varp] = list()
    RMSD['density'] = list()
    for varc in ['nittemp','nitdens']:
        corrfloat[varc] = list()
        corrmod[varc] = list()
        covnfloat[varc] = list()
        covnmod[varc] = list()
    for ip, p in enumerate(list_float_track):
        Lon[ip] = p.lon
        Lat[ip] = p.lat

#        timelabel_list.append(p.time)
#        for var_mod in ['N3n','Chla']:
#        for var_mod in ['N3n','P_l']:
        for var_mod in ['N3n']:
            if p.available_params.find(FLOATVARS[var_mod])<0 : continue
            timelabel_list.append(p.time)
            NewProf_5m = np.zeros(len(NewPres_5m))
            NewProf_5m[:] = np.nan
            # print p.time
            # print var_mod
            var = FLOATVARS[var_mod]
            adj=Adj[var_mod]

            # FLOAT
            Pres,Prof,Qc = p.read(var,read_adjusted=adj)
            ii = Pres<=300
            # ii = Pres<=200
            # Deve avere almeno 5 records:
            if len(Prof[ii])>3 :
                NewProf_5m = np.interp(NewPres_5m,Pres[ii],Prof[ii])

            # MODEL
            TM = MM.modeltime(p)
            FILENAME = BASEDIR + TM.strftime("PROFILES/ave.%Y%m%d-%H:00:00.profiles.nc")
            M = readModelProfile(FILENAME,var_mod,p.ID())
            M_newDepth=np.interp(NewPres_5m,TheMask.zlevels[:max_depth+1],M[:max_depth+1])

            rmsdfloat = (np.sum((M_newDepth-NewProf_5m)**2)/nPnewpres)**.5
            RMSD[var_mod].append(rmsdfloat)

            MP_newDepth = {}
            NewProfp_5m = {}
            if var_mod == 'N3n':
                for ivarp, varp_mod in enumerate(VARLIST_PHYS):
        #     if (ivar==0):
                    varp = FLOATVARS[varp_mod]
                    adj=Adj[varp_mod]

                    Presp,Profp,Qcp = p.read(varp,read_adjusted=adj)
                    ii = Presp<=300
                    # Deve avere almeno 5 records:
                    if len(Profp[ii])>5 :
                        NewProfp_5m[varp_mod] = np.interp(NewPres_5m,Presp[ii],Profp[ii])


        # PLOT FOR THE MODEL
                    TMP = MM.modeltime(p)
                    FILENAME = BASEDIR + TMP.strftime("PROFILES/ave.%Y%m%d-%H:00:00.profiles.nc")
                    MP = readModelProfile(FILENAME,varp_mod,p.ID())
                    MP_newDepth[varp_mod]=np.interp(NewPres_5m,TheMask_phys.zlevels[:max_depthp+1],MP[:max_depthp+1])

                    rmsdfloat = (np.sum((MP_newDepth[varp_mod]-NewProfp_5m[varp_mod])**2)/nPnewpres)**.5
                    RMSD[varp_mod].append(rmsdfloat)

                MP_newDepth['density'] = seawater.dens(MP_newDepth['vosaline'],MP_newDepth['votemper'],NewPres_5m)
                NewProfp_5m['density'] = seawater.dens(NewProfp_5m['vosaline'],NewProfp_5m['votemper'],NewPres_5m)

                rmsdfloat = (np.sum((MP_newDepth['density']-NewProfp_5m['density'])**2)/nPnewpres)**.5
                RMSD['density'].append(rmsdfloat)

                corrcoef = np.corrcoef(NewProfp_5m['votemper'],NewProf_5m)
                corrfloat['nittemp'].append(corrcoef[0,1])
                corrcoef = np.corrcoef(NewProfp_5m['density'],NewProf_5m)
                corrfloat['nitdens'].append(corrcoef[0,1])

                covmat = np.cov(NewProfp_5m['votemper'],NewProf_5m)
                covnfloat['nittemp'].append(covmat[0,1]/covmat[0,0])
                covmat = np.cov(NewProfp_5m['density'],NewProf_5m)
                covnfloat['nitdens'].append(covmat[0,1]/covmat[0,0])

                corrcoef = np.corrcoef(MP_newDepth['votemper'],M_newDepth)
                corrmod['nittemp'].append(corrcoef[0,1])
                corrcoef = np.corrcoef(MP_newDepth['density'],M_newDepth)
                corrmod['nitdens'].append(corrcoef[0,1])

                covmat = np.cov(MP_newDepth['votemper'],M_newDepth)
                covnmod['nittemp'].append(covmat[0,1]/covmat[0,0])
                covmat = np.cov(MP_newDepth['density'],M_newDepth)
                covnmod['nitdens'].append(covmat[0,1]/covmat[0,0])


    LIST = [ii for ii in range(4)]
    LIST[0] = timelabel_list
    LIST[1] = corrmod
    LIST[2] = covnmod
    LIST[3] = RMSD

    if p.available_params.find(FLOATVARS[var_mod])>0:
      filename = OUTDIR + 'corrcov' + run + '_' + p.name() + '.pkl'
      fid = open(filename,'wb')
      pickle.dump(LIST,fid)
      fid.close()

      if dofloat==True:
        LISTF = [ii for ii in range(3)]
        LISTF[0] = timelabel_list
        LISTF[1] = corrfloat
        LISTF[2] = covnfloat

        filename = OUTDIR + 'corrcovFloat_' + p.name() + '.pkl'
        fid = open(filename,'wb')
        pickle.dump(LISTF,fid)
        fid.close()



import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Exclusion of profiles not consistent with ref at depth
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    # parser.add_argument(   '--inwoa', '-w',
    #                             type = str,
    #                             default = None,
    #                             required = True,
    #                             help = "Directory with files of WOA")

    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    parser.add_argument(   '--topdepth','-t',
                                type = str,
                                required = True,
                                help = 'Top of layer on which evaluate RMSD')
    parser.add_argument(   '--bottomdepth','-b',
                                type = str,
                                required = True,
                                help = 'Bottom of layer on which evaluate RMSD')   
    return parser.parse_args()

args = argument()



import numpy as np
import pylab as pl
import basins.OGS as basV2
import figure_generator
import netCDF4
import scipy.io.netcdf as NC
import datetime
from instruments import lovbio_float as bio_float
from static.climatology import get_climatology
from basins.basin import ComposedBasin
# from common.dataextractor import DataExtractor
from commons.layer import Layer
from commons.mask import Mask
from commons.submask import SubMask
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.utils import addsep
from instruments.matchup_manager import Matchup_Manager
from instruments.var_conversions import LOVFLOATVARS
from profilerDAf import BASEDIR,TL,ALL_PROFILES,T_INT
from timeseries.plot import Hovmoeller_matrix
from timeseries.plot import read_pickle_file, read_basic_info


IDrun='floatcfr'
OUTDIR=addsep(args.outdir)
# MODDIR=addsep(args.inmodel)

TheMask = Mask(args.maskfile)

def readModelProfile(filename,var, wmo):
    ncIN = NC.netcdf_file(filename,'r')
    M = ncIN.variables[var].data.copy()
    iProfile = ncIN.CruiseIndex.rsplit(", ").index(wmo)
    ncIN.close()
    Profile = M[iProfile,:]
    return Profile


# TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")
# maskfile8="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/V1/meshmask_872.nc"
# Mask8 = Mask(maskfile8)
# jpk8,jpj8,jpi8 = Mask8.shape
# TheMask= Mask(args.maskfile, loadtmask=False)
# jpk,jpj,jpi = TheMask.shape
# z = -TheMask.zlevels

LIMdep = [np.float(args.topdepth),np.float(args.bottomdepth)]
txtLIMdep = '%d_' %LIMdep[0] + '%d' %LIMdep[1]

PresDOWN=np.array([0,25,50,75,100,125,150,200,400,600,800,1000])
LayerList=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]

NewPres_5m=np.linspace(0,40,9)
NewPres_10m=np.linspace(50,200,(200-50)/10+1)
NewPres_25m=np.linspace(225,1000,(1000-225)/25+1)

NewPres = np.concatenate([NewPres_5m,NewPres_10m,NewPres_25m])
maskPres = (NewPres>LIMdep[0]) & (NewPres<=LIMdep[1])

# z_clim = np.array([(l.bottom+l.top)/2  for l in LayerList])
max_depth = TheMask.getDepthIndex(PresDOWN[-1])+1


# TL = TimeList.fromfilenames(TI, MODDIR, "ave*nc")
Profilelist_1=bio_float.FloatSelector(None,T_INT,basV2.med)
wmo_list=bio_float.get_wmo_list(Profilelist_1)

# SUBLIST = basV2.P.basin_list
# nSub = len(SUBLIST)


# N3n_clim, N3n_std = get_climatology('N3n', SUBLIST, LayerList)
# N1p_clim, N1p_std = get_climatology('N1p', SUBLIST, LayerList)
# O2o_clim, O2o_std = get_climatology('O2o', SUBLIST, LayerList)



# VARLIST=['P_l','N1p','N3n','O2o']
VARLIST=['P_l','N3n']
# var_dype = [(var,np.float32) for var in VARLIST]
# nVar = len(VARLIST)
limVAR = {
    'P_l': 2,
    'N3n': .75,
}

Adj = {
    'P_l': True,
    'N3n': True,
    'O2o': False,
}

nZ = NewPres.shape[0]
zFloat = -NewPres

DICTiSub = {}
for iSub,sub in enumerate(basV2.P):
    DICTiSub[sub.name] = iSub

PresVar = {}
for var_mod in VARLIST:
    PresVar[var_mod] = []

# LIST = {}
Nover = {}
for var_mod in VARLIST:
    Nover[var_mod] = 0
    # LIST[var_mod] = [[] for ii in range(3)]

thresholdRMSD = 1.5


MM = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

# for wmo in [wmo_list[0]]:
RMSDmasked = {}

# deptlim = 0
for wmo in wmo_list:
    RMSDmasked[wmo] = []
    # iip = {}
    # iipcheck = {}
    # for var_mod in VARLIST:
    #     iip[var_mod] = 0
    #     iipcheck[var_mod] = 0
    print wmo
    list_float_track=bio_float.filter_by_wmo(Profilelist_1,wmo)
    np_tot = len(list_float_track)
    print '  number of profiles for this float ' + np.str(np_tot) 
        
    for p in list_float_track:
        # for indSub, sub in enumerate(basV2.Pred):
        timef12 = datetime.datetime(p.time.year,p.time.month,p.time.day,12,0)
        # for indSub,sub in enumerate(basV2.Pred):
        #     if sub.is_inside(p.lon,p.lat):
        #         iSubp = indSub
        if(not(timef12 in TL.Timelist)): continue
        for var_mod in VARLIST:
            var = LOVFLOATVARS[var_mod]
            adj=Adj[var_mod]
            Pres,Prof,Qc=p.read(var,read_adjusted=adj)
            # if Pres.size>0:
            #     PresVar[var_mod].append(Pres)
            ii = Pres<=1000 # Deve avere almeno 5 records:
            if len(Prof[ii])>5:
                if (Prof[0])>limVAR[var_mod]: Nover[var_mod]+=1
                # iip[var_mod] += 1
                if var_mod=='P_l': continue
                NewProf = np.interp(NewPres,Pres[ii],Prof[ii])

                TM = MM.modeltime(p)
                FILENAME = BASEDIR + TM.strftime("PROFILES/ave.%Y%m%d-%H:00:00.profiles.nc")
                M = readModelProfile(FILENAME,var_mod,p.ID())
                NewModel = np.interp(NewPres,TheMask.zlevels[:max_depth],M[:max_depth])

                RMSDmasked[wmo].append((np.nanmean((NewProf[maskPres]-NewModel[maskPres])**2))**.5)
                # NewClim = np.interp(NewPres,z_clim,N3n_clim[iSubp,:])
                # if RMSDmasked>thresholdRMSD :
                #     # print '    ...  excluding a profile for comparison with clim'
                #     iipcheck[var_mod] += 1



pl.close('all')

prctileLIST = [95,75,50,25,5]
DICTprctile = {}
for prc in prctileLIST:
    DICTprctile[prc] = []

wmo_plotlist = []
for wmo in wmo_list:
    if RMSDmasked[wmo]==[]: continue
    wmo_plotlist.append(wmo)
    for prc in prctileLIST:
        DICTprctile[prc].append(np.nanpercentile(RMSDmasked[wmo],prc))


pl.figure(figsize=(9,5))

for prc in prctileLIST:
    pl.bar(wmo_plotlist,DICTprctile[prc],1,label=prc)

pl.ylim(0,2.5)
pl.tick_params(axis='x',rotation=30)
ax = pl.gca()
for ll in ax.xaxis.get_ticklabels():
    ll.set_horizontalalignment('right')
pl.title(txtLIMdep)

pl.legend()
pl.grid()

pl.savefig(OUTDIR + '/percentileRMSD' + txtLIMdep + '.png')


thresholdLIST = [.5,1,1.5,2]
DICTthresh = {}
for trs in thresholdLIST:
    DICTthresh[trs] = []

wmo_plotlist = []
wmo_totp = []
for wmo in wmo_list:
    if RMSDmasked[wmo]==[]: continue
    wmo_plotlist.append(wmo)
    wmo_totp.append(len(RMSDmasked[wmo]))
    for trs in thresholdLIST:
        DICTthresh[trs].append(np.nansum(np.array(RMSDmasked[wmo])>trs))

nTOTp = np.sum(wmo_totp)

pl.figure(figsize=(9,5))

for trs in thresholdLIST:
    pl.bar(wmo_plotlist,DICTthresh[trs],1, \
        label=np.str(trs) + \
        ' ' + np.str(np.sum(DICTthresh[trs])) + ' on ' + np.str(nTOTp) \
        )

pl.tick_params(axis='x',rotation=30)
ax = pl.gca()
for ll in ax.xaxis.get_ticklabels():
    ll.set_horizontalalignment('right')
pl.title(txtLIMdep)

pl.legend()
pl.grid()

pl.savefig(OUTDIR + '/thresholdRMSD' + txtLIMdep + '.png')


pl.figure(figsize=(9,5))

for trs in thresholdLIST:
    pl.bar(wmo_plotlist,100.*np.array(DICTthresh[trs])/np.array(wmo_totp),1, \
        label=np.str(trs) + \
        ' ' + np.str(np.sum(DICTthresh[trs])) + ' on ' + np.str(nTOTp) \
        )

pl.tick_params(axis='x',rotation=30)
ax = pl.gca()
for ll in ax.xaxis.get_ticklabels():
    ll.set_horizontalalignment('right')
pl.title(txtLIMdep)

pl.legend()
pl.grid()

pl.savefig(OUTDIR + '/thresholdpercRMSD' + txtLIMdep + '.png')


pl.show(block=False)



#     for var_mod in VARLIST:
#         if iip[var_mod]==0:
#             print '    Any profile for ' + var_mod
#         else:
#             LIST[var_mod][0].append(wmo[-3:])
#             LIST[var_mod][1].append(iip[var_mod])
#             LIST[var_mod][2].append(iip[var_mod]-iipcheck[var_mod])

#             percex = np.float(iipcheck[var_mod])/np.float(iip[var_mod])*100.
#             print '  Number of profiles for ' + var_mod + \
#                 ' ' + np.str(iip[var_mod])
#             print '  Number of excluded ' + \
#                 '  ' + np.str(iipcheck[var_mod]) + ' percent ' + \
#                 np.str(percex) + '  (considered ' + np.str(100.-percex) + ')'


# pl.close('all')
# fig,axs = pl.subplots(3,1,figsize=[12,8],sharex=True)

# pl.sca(axs[0])
# pl.title('Number of profiles available for N3n')
# pl.bar(LIST['N3n'][0],LIST['N3n'][1])
# pl.grid()

# pl.sca(axs[1])
# pl.title('Number of non-excluded profiles')
# pl.bar(LIST['N3n'][0],LIST['N3n'][2])
# pl.grid()

# pl.sca(axs[2])
# pl.title('Percentage of non-excluded profiles')
# percnoex = np.array(LIST['N3n'][2]).astype(float) / np.array(LIST['N3n'][1]).astype(float) * 100.
# pl.bar(LIST['N3n'][0],percnoex)
# pl.grid()

# pl.show(block=False)
# pl.savefig('excluded_profiles_2015.png')
    
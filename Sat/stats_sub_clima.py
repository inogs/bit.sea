import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Produces statistics (mean, std and others) over Mediterranean subbasins
    from a climatology
    (to be launche once for each climatology)
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(   '--climfile', '-c',
                                type = str,
                                required = True,
                                help = ''' Climatology .nc file used to apply check on sat data, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/CCI02/SatClimatology.nc'''

                                )
    
    parser.add_argument(   '--mesh', '-m',
                                type = str,
                                required = True,
                                choices = ['SatOrigMesh','V4mesh','V1mesh','KD490mesh','SAT1km_mesh', 'Mesh24'],
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
                                )

    parser.add_argument(   '--submaskdir', '-s',
                                type = str,
                                required = True,
                                help = '''  Directory with subbasins on satellite mask'''
                                )

    parser.add_argument(   '--statsdir', '-w',
                                type = str,
                                required = True,
                                help = ''' STATS  directory'''

                                )

    return parser.parse_args()


args = argument()
import pickle
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins import V2
from postproc import masks
import numpy as np
from commons.utils import addsep
import os

import SatManager as Sat


SUBMASKDIR  = addsep(args.submaskdir)
CLIM_FILE = args.climfile
STATSDIR  = addsep(args.statsdir)

maskSat = getattr(masks,args.mesh)

reset = False

somestats = False
for ii in range(365):
    strday = "%03d" %(ii+1)
    climstatsfile = STATSDIR + 'statsday' + strday + '_clim.pkl'
    if (not os.path.exists(climstatsfile)):
        somestats = True
        break

if somestats or reset:
    print('Read climatology')
    MEAN,STD = Sat.readClimatology(CLIM_FILE)

    filemaskmed = SUBMASKDIR + 'maskmed_1km.npy'
    maskmed_1km = np.load(filemaskmed)

    masksub_M = {}
    for sub in V2.P:
        filemasksub = SUBMASKDIR + 'masksub.' + sub.name + '.npy'
        masksub = np.load(filemasksub)
        masksub_M[sub.name] = np.zeros((maskSat.jpj,maskSat.jpi),dtype=bool)
        masksub_M[sub.name][maskmed_1km] = masksub


else:
    print('All climatology statistics done')


for ii in range(365):
    strday = "%03d" %(ii+1)
    filestatsclim = STATSDIR + 'statsday' + strday + '_clim.pkl'
    exit_condition = os.path.exists(filestatsclim)

    if (exit_condition) and (reset==False):
        continue

    julian = ii+1
    print(' ... day ' + np.str(julian))

    DAILY_REF_MEAN = MEAN[julian-1,:,:]
    DAILY_REF_STD  =  STD[julian-1,:,:]    
    
    stats_clima = {}
    for sub in V2.P:
        stats_clima[sub.name] = np.zeros(12)
        stats_clima[sub.name][:] = np.nan

        climadmean = DAILY_REF_MEAN[masksub_M[sub.name]]
        climadmean[climadmean<0] = np.nan
        climadstd = DAILY_REF_STD[masksub_M[sub.name]]
        climadstd[climadmean<0] = np.nan

        climadmadd_2std = climadmean + 2*climadstd
        climadmsub_2std = climadmean - 2*climadstd


        if np.nansum(masksub_M[sub.name])>0:
            stats_clima[sub.name][0] = np.nanmean(climadmean)
            stats_clima[sub.name][1] = np.nanmin(climadmean)
            stats_clima[sub.name][2] = np.nanmax(climadmean)
            stats_clima[sub.name][3] = np.nanmean(climadmadd_2std)
            stats_clima[sub.name][4] = np.nanmin(climadmadd_2std)
            stats_clima[sub.name][5] = np.nanmax(climadmadd_2std)
            stats_clima[sub.name][6] = np.nanmean(climadmsub_2std)
            stats_clima[sub.name][7] = np.nanmin(climadmsub_2std)
            stats_clima[sub.name][8] = np.nanmax(climadmsub_2std)
            stats_clima[sub.name][9] = np.nanmean(climadstd)
            stats_clima[sub.name][10] = np.nanmin(climadstd)
            stats_clima[sub.name][11] = np.nanmax(climadstd)



    fid = open(filestatsclim,'wb')
    pickle.dump(stats_clima,fid)
    fid.close()

    print 'Done statistics with file ' + strday + ' of 365 ' 
    print '   ---------------------------------------------------   '


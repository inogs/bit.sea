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

    parser.add_argument(   '--masktype', '-t',
                                type = str,
                                required = True,
                                choices = ['Open','All'],
                                help = ''' Type of mask (Open or All)'''

                                )

    return parser.parse_args()


args = argument()
import scipy.io.netcdf as NC
from basins import V2
from postproc import masks
import numpy as np
from commons.utils import addsep
import os

import SatManager as Sat


SUBMASKDIR  = addsep(args.submaskdir)
CLIM_FILE = args.climfile
STATSDIR  = addsep(args.statsdir)
masktype = args.masktype

maskSat = getattr(masks,args.mesh)

reset = True

somestats = False
for ii in range(365):
    strday = "%03d" %(ii+1)
    climstatsfile = STATSDIR + 'statsday' + strday + '_clim' + masktype + '.nc'
    if (not os.path.exists(climstatsfile)):
        somestats = True
        break

if somestats or reset:
    print('Read climatology and mask for subbasins')
    MEAN,STD = Sat.readClimatology(CLIM_FILE)

    filemaskmed = SUBMASKDIR + 'maskmed_1km' + masktype + '.npy'
    maskmed_1km = np.load(filemaskmed)

    masksub_M = {}
    nsub = 0
    subnames = ''
    for sub in V2.P:
        print sub.name
        nsub += 1
        subnames += sub.name + ', '
        filemasksub = SUBMASKDIR + 'masksub.' + sub.name + masktype + '.npy'
        masksub = np.load(filemasksub)
        masksub_M[sub.name] = np.zeros((maskSat.jpj,maskSat.jpi),dtype=bool)
        masksub_M[sub.name][maskmed_1km] = masksub




else:
    print('All climatology statistics done')


for ii in range(365):
    strday = "%03d" %(ii+1)
    filestatsclim = STATSDIR + 'statsday' + strday + '_clim' + masktype + '.nc'
    exit_condition = os.path.exists(filestatsclim)

    if (exit_condition) and (reset==False):
        continue

    julian = ii+1
    print(' ... day ' + np.str(julian))

    DAILY_REF_MEAN = MEAN[julian-1,:,:]
    DAILY_REF_STD  =  STD[julian-1,:,:]    
    
    stats_clima = np.zeros((nsub,13))
    stats_clima[:,:] = np.nan
    for isub,sub in enumerate(V2.P):

        climadmean = DAILY_REF_MEAN[masksub_M[sub.name]]
        climadmean[climadmean<0] = np.nan
        climadstd = DAILY_REF_STD[masksub_M[sub.name]]
        climadstd[np.isnan(climadmean)] = np.nan

        climadmadd_2std = climadmean + 2*climadstd
        climadmsub_2std = climadmean - 2*climadstd


        if np.nansum(masksub_M[sub.name])>0:
            stats_clima[isub,0] = np.nanmean(climadmean)
            stats_clima[isub,1] = np.nanmin(climadmean)
            stats_clima[isub,2] = np.nanmax(climadmean)
            stats_clima[isub,3] = np.nanmean(climadmadd_2std)
            stats_clima[isub,4] = np.nanmin(climadmadd_2std)
            stats_clima[isub,5] = np.nanmax(climadmadd_2std)
            stats_clima[isub,6] = np.nanmean(climadmsub_2std)
            stats_clima[isub,7] = np.nanmin(climadmsub_2std)
            stats_clima[isub,8] = np.nanmax(climadmsub_2std)
            stats_clima[isub,9] = np.nanmean(climadstd)
            stats_clima[isub,10] = np.nanmin(climadstd)
            stats_clima[isub,11] = np.nanmax(climadstd)
            stats_clima[isub,12] = np.nanstd(climadmean)


    ncOUT = NC.netcdf_file(filestatsclim,'w')
    ncOUT.createDimension('subbasin',nsub)
    ncOUT.createDimension('stattype',13)

    setattr(ncOUT,'sublist',subnames)

    ncvar = ncOUT.createVariable('SUBstatistics','f',('subbasin','stattype'))
    ncvar[:] = stats_clima

    ncOUT.close()


    print 'Done statistics with file ' + strday + ' of 365 ' 
    print '   ---------------------------------------------------   '


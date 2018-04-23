import numpy as np
import argparse
from Sat import SatManager as Sat
#from commons.time_interval import TimeInterval
#from commons.Timelist import TimeList
from basins import V2
from postproc import masks
from commons.utils import addsep


def argument():
    parser = argparse.ArgumentParser(description = '''
    Apply check based on climatology to sat ORIG files
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' CHECKED sat directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/'''

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

    return parser.parse_args()

args = argument()


OUTDIR  = addsep(args.outdir)
CLIM_FILE = args.climfile

maskSat = getattr(masks,args.mesh)


readClimatology = True
if readClimatology:
    print('--- Read climatology')
    MEAN,STD = Sat.readClimatology(CLIM_FILE)
    # Build a Mediterranean mask
    mask365 = MEAN>0
    maskmed_1km = np.sum(mask365,0)>0

filename = OUTDIR + 'maskmed_1km.npy'
print(filename)
np.save(filename,maskmed_1km)

# Extracting points with observations
indlon = range(maskSat.jpi)
indlat = range(maskSat.jpj)
indlonM,indlatM = np.meshgrid(indlon,indlat)
indlon_masked = indlonM[maskmed_1km]
indlat_masked = indlatM[maskmed_1km]

LonM,LatM = np.meshgrid(maskSat.lon,maskSat.lat)
Lon_masked = LonM[maskmed_1km]
Lat_masked = LatM[maskmed_1km]

numP = Lat_masked.shape[0]

print('   Cycle on sub')
for sub in V2.Pred:
    masksub = np.zeros(numP,dtype=bool)
    print(sub.name)
    for iip in range(numP):
        lonp = Lon_masked[iip]
        latp = Lat_masked[iip]
        indlonp = indlon_masked[iip]
        indlatp = indlat_masked[iip]
        if sub.is_inside(lonp,latp):
            masksub[iip] = True

    Np_sub = np.sum(masksub)
    print('     ... Np ' + np.str(Np_sub))

    filename = OUTDIR + 'masksub.' + sub.name + '.npy'
    print(filename)
    np.save(filename,masksub)


masksub = np.zeros(numP,dtype=bool)
masksub[:] = True
Np_sub = np.sum(masksub)
print('med')
print('     ... Np ' + np.str(Np_sub))

filename = OUTDIR + 'masksub.med.npy'
print(filename)
np.save(filename,masksub)






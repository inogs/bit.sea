import numpy as np
import argparse
from Sat import SatManager as Sat
#from commons.time_interval import TimeInterval
#from commons.Timelist import TimeList
from basins import V2
from postproc import masks
from commons.utils import addsep
from commons.mask import Mask
from Sat import interp2d


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

    parser.add_argument(   '--modelmeshfile', '-l',
                                type = str,
                                required = True,
                                help = ''' File of the model mesh '''
                                )

    return parser.parse_args()

args = argument()


OUTDIR  = addsep(args.outdir)
CLIM_FILE = args.climfile
Maskfile = args.modelmeshfile

maskSat = getattr(masks,args.mesh)
maskMod = Mask(Maskfile)


readClimatology = True
readClimatology = False
if readClimatology:
    print '--- Read climatology'
    MEAN,STD = Sat.readClimatology(CLIM_FILE)
    # Build a Mediterranean mask
    mask365 = np.zeros_like(MEAN)
    mask365[MEAN>0] = 1
    maskmed_1km = (np.sum(mask365,0))>0

    # Indexes for interpolation
    print '--- Indexes for interpolation '
    x = maskMod.xlevels[0,:]
    y = maskMod.ylevels[:,0]
    x1km = maskSat.lon
    y1km = maskSat.lat
    I_START, I_END = interp2d.array_of_indices_for_slicing(x, x1km)
    J_START, J_END = interp2d.array_of_indices_for_slicing(y, y1km)


# Considering coastal areas
maskmodel200 = maskMod.mask_at_level(200)
maskmodel0 = maskMod.mask_at_level(0)
maskmodelcoast = maskmodel0.copy()
maskmodelcoast[maskmodel200] = False

indlonCoast,indlatCoast = np.nonzero(maskmodelcoast)
Ncoast = indlonCoast.shape[0]
for iip in range(Ncoast):
    maskmed_1km[J_START[indlonCoast[iip]]:J_END[indlonCoast[iip]],
                I_START[indlatCoast[iip]]:I_END[indlatCoast[iip]]] = False
    


filename = OUTDIR + 'maskmed_1kmOpen.npy'
print(filename)
np.save(filename,maskmed_1km)

# Extracting points with observations and open sea
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
    if sub.name in ['lev1','lev2','lev3','lev4']:
        masksub = np.zeros(numP,dtype=bool)
        print(sub.name)
        for iip in range(numP):
            lonp = Lon_masked[iip]
            latp = Lat_masked[iip]
            # indlonp = indlon_masked[iip]
            # indlatp = indlat_masked[iip]
            if sub.is_inside(lonp,latp):
                masksub[iip] = True

        Np_sub = np.sum(masksub)
        print('     ... Np ' + np.str(Np_sub))

        filename = OUTDIR + 'masksub.' + sub.name + 'Open.npy'
        print(filename)
        np.save(filename,masksub)


masksub = np.zeros(numP,dtype=bool)
masksub[:] = True
Np_sub = np.sum(masksub)
print('med')
print('     ... Np ' + np.str(Np_sub))

filename = OUTDIR + 'masksub.medOpen.npy'
print(filename)
np.save(filename,masksub)






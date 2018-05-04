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

    parser.add_argument(   '--masktype', '-t',
                                type = str,
                                required = True,
                                choices = ['Open','All'],
                                help = ''' If Open or All mask '''
                                )

    return parser.parse_args()

args = argument()


OUTDIR  = addsep(args.outdir)
CLIM_FILE = args.climfile
masktype = args.masktype
Maskfile = args.modelmeshfile

maskSat = getattr(masks,args.mesh)
maskMod = Mask(Maskfile)


readClimatology = True
if readClimatology:
    print '--- Read climatology'
    MEAN,STD = Sat.readClimatology(CLIM_FILE)
    # Build a Mediterranean mask
    mask365 = np.zeros_like(MEAN)
    mask365[MEAN>0] = 1
    maskmed_1km = (np.sum(mask365,0))>0
    if masktype=='Open':
        masksatOpen = np.zeros_like(maskmed_1km)

    # Indexes for interpolation
    print '--- Indexes for interpolation '
    x = maskMod.xlevels[0,:]
    y = maskMod.ylevels[:,0]
    x1km = maskSat.lon
    y1km = maskSat.lat
    I_START, I_END = interp2d.array_of_indices_for_slicing(x, x1km)
    J_START, J_END = interp2d.array_of_indices_for_slicing(y, y1km)


# Considering coastal areas

if masktype=='Open':
    maskmodel200 = maskMod.mask_at_level(200)
    indlonOpen,indlatOpen = np.nonzero(maskmodel200)
    Nopen = indlonOpen.shape[0]
    for iip in range(Nopen):
        masksatOpen[J_START[indlonOpen[iip]]:J_END[indlonOpen[iip]],
                    I_START[indlatOpen[iip]]:I_END[indlatOpen[iip]]] = True
    
    maskmed_1km[masksatOpen==False] = False


filename = OUTDIR + 'maskmed_1km' + masktype + '.npy'
print(filename)
np.save(filename,maskmed_1km)

# Extracting points with observations and in the mask (open or all)
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
        if sub.is_inside(lonp,latp):
            masksub[iip] = True

    Np_sub = np.sum(masksub)
    print('     ... Np ' + np.str(Np_sub))

    filename = OUTDIR + 'masksub.' + sub.name + masktype + '.npy'
    print(filename)
    np.save(filename,masksub)


masksub = np.zeros(numP,dtype=bool)
masksub[:] = True
Np_sub = np.sum(masksub)
print('med')
print('     ... Np ' + np.str(Np_sub))

filename = OUTDIR + 'masksub.med' + masktype + '.npy'
print(filename)
np.save(filename,masksub)






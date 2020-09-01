import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Calculate var_mod as difference between vardiff and varsat
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Outputdir'''
                            )

    parser.add_argument(   '--inmod', '-d',
                            type = str,
                            required = True,
                            help = 'Input model'
                            )

    parser.add_argument(   '--insat', '-s',
                            type = str,
                            required = True,
                            help = 'Input var sat'
                            )

#    parser.add_argument(   '--limgib', '-g',
#                            type = str,
#                            required = True,
#                            help = 'Liongitude limit of Gibraltar'
#                            )

    parser.add_argument(   '--maskfile', '-m',
                            type = str,
                            required = True,
                            help = 'model maskfile meshmask.nc'
                            )

    parser.add_argument(   '--minvmod', '-v',
                            type = str,
                            required = True,
                            nargs=2,
                            help = 'min mod var quotas of vdiff in open sea and coast'
                            )

    return parser.parse_args()


args = argument()

try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False


import numpy as np
from commons.dataextractor import DataExtractor
from commons.mask import Mask 
from commons.Timelist import TimeList
from commons.utils import addsep
from commons.timerequestors import Clim_month
from postproc import masks
from Sat import SatManager

maskSat = getattr(masks,'Mesh24')

INSAT = addsep(args.insat)
INMOD = addsep(args.inmod)
OUTDIR = addsep(args.outdir)
TheMask = Mask(args.maskfile)
minvmodLIST = args.minvmod

_,jpj,jpi = TheMask.shape
indGib,_ = TheMask.convert_lon_lat_to_indices(-5.2,35.9)

masksurf = TheMask.mask[0,:,:].copy()
masksurf[:,:indGib] = False

maskOpen = (TheMask.mask_at_level(200)) & (masksurf)
maskCoast =(TheMask.mask_at_level(200)==False) & (masksurf) 

minOpen = float(minvmodLIST[0])
minCoast = float(minvmodLIST[1])

TL_mod = TimeList.fromfilenames(None,INMOD,'ave*Chla.nc',prefix='ave.')


MONTHLY_reqs=[Clim_month(i) for i in range(1,13)]
Chl = np.zeros((12,jpj,jpi), dtype=float)
ChlSquare = np.zeros((12,jpj,jpi), dtype=float)

for req in MONTHLY_reqs[rank::nranks]:
    im = req.month
    print 'Month %i' %im
    filesat = INSAT + "/var2Dsat.CCI.F7.2.%02d.nc"%(im)
    De = DataExtractor(TheMask,filename=filesat,varname='variance',dimvar=2)
    varSat = De.filled_values[:,:]

    ii,w = TL_mod.select(req)
    nFiles = len(ii)
    M  = np.zeros((nFiles,jpj,jpi),np.float32)
    M2 = np.zeros((nFiles,jpj,jpi),np.float32)

    for iFrame, j in enumerate(ii):
        inputfile = TL_mod.filelist[j]
        print '   ' + inputfile
        De = DataExtractor(TheMask,filename=inputfile,varname='Chla',dimvar=3)
        CHL = De.filled_values[0,:,:].copy()
        TmpMat = np.zeros(CHL.shape)
        M[iFrame,:,:]  = CHL

        TmpMat = CHL*CHL
        M2[iFrame,:,:] = TmpMat

    MonthIndex = req.month-1
    Chl[MonthIndex,:,:] = np.nanmean(M,0)
    ChlSquare[MonthIndex,:,:] = np.nanmean(M2,0)

    varMod = ChlSquare[MonthIndex,:,:] - Chl[MonthIndex,:,:]*Chl[MonthIndex,:,:]


    MonthStr = '%02d'%(im)

    ratio = varMod/varSat
    fname = OUTDIR + '/ratiovarM_S.' + MonthStr + '.nc'
    print "\tsaving ", fname
    SatManager.dumpGenericNativefile(fname,ratio,'variance',maskSat)

    fname = OUTDIR + '/var2D.' + MonthStr + '.nc'
    print "\tsaving ", fname
    SatManager.dumpGenericNativefile(fname,varMod,'variance',maskSat)

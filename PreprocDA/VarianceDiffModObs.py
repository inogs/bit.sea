import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Computes the monthly variance of the difference between Model output and Satellite observations
  
    MyMesh must be consistent with both the Model data inputs (within ModDir)
    and the Satellite observations  (within SatDir)

    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Outputdir'''
                            )

    parser.add_argument(   '--insat', '-s',
                            type = str,
                            required = True,
                            help = 'Input map sat'
                            )

    parser.add_argument(   '--inmod', '-d',
                            type = str,
                            required = True,
                            help = 'Input mod chl'
                            )

    parser.add_argument(   '--mesh', '-m',
                            type = str,
                            required = True,
                            help = 'Mesh of sat map resolution'
                            )

    parser.add_argument(   '--maskfile', '-f',
                            type = str,
                            required = True,
                            help = 'model maskfile meshmask.nc'
                            )

    return parser.parse_args()

args = argument()


import numpy as np
from commons.dataextractor import DataExtractor
from commons.mask import Mask
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.timerequestors import Clim_month
from commons.utils import addsep
from postproc import masks
import Sat.SatManager as Sat
import netCDF4

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

ModDir = addsep(args.inmod)
SatDir = addsep(args.insat)
OUTDIR = addsep(args.outdir)
MyMesh = getattr(masks,args.mesh)
TheMask = Mask(args.maskfile)


TLModel = TimeList.fromfilenames(None, ModDir,"*Chla.nc",prefix='ave.')
Timestart = TLModel.Timelist[0].strftime('%Y%m%d')
Time__end = TLModel.Timelist[-1].strftime('%Y%m%d')
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TLSat   = TimeList.fromfilenames(TI, SatDir,'*.nc',prefix='',dateformat='%Y%m%d')
MONTHLY_reqs=[Clim_month(i) for i in range(1,13)]

jpi = MyMesh.jpi
jpj = MyMesh.jpj

ChlDiffSquared = np.zeros((12,jpj,jpi), dtype=float)

# get fillValue assuming that all files in the
# folder have the same fillValue
ncIN = netCDF4.Dataset(TLSat.filelist[0],'r')
varObj = ncIN.variables['CHL']
fillValue = varObj.getncattr('_FillValue')
ncIN.close()


for req in MONTHLY_reqs[rank::nranks]:
    print req
    iiMod, wMod = TLModel.select(req)
    #iiSat, wSat = TLSat.select(req)

    #assert(len(iiMod) == len(iiSat)), "TimeLists not coherent"

    nFiles = len(iiMod)
    DiffModSat = np.zeros((nFiles,jpj,jpi))
    DiffModSat[:,:,:] = np.nan

    for iFrame, j in enumerate(iiMod):

        ModFile = TLModel.filelist[j]
        ModDate = TLModel.Timelist[j]
        if (TLSat.Timelist[0]-ModDate).days>3: continue
        if (ModDate-TLSat.Timelist[-1]).days>3: continue
        print ModDate
        jsat = TLSat.find(ModDate)
        SatFile = TLSat.filelist[jsat]

        SatChl = Sat.readfromfile(SatFile,'CHL')
        De = DataExtractor(TheMask,filename=ModFile,varname='Chla',dimvar=3)
        ModChl = De.filled_values[0,:,:]


        MyMask = (SatChl == fillValue) | (np.isnan(ModChl))

        # debug output
        # print "TLModel.Timelist[j] ", TLModel.Timelist[j], "  TLSat.Timelist[j] ",TLSat.Timelist[j], " j = ", j

        assert(DiffModSat[0,:,:].shape == ModChl.shape)
        assert(DiffModSat[0,:,:].shape == SatChl.shape)

        TmpMatrix = ModChl - SatChl
        TmpMatrix = TmpMatrix*TmpMatrix      
        TmpMatrix[MyMask] = np.nan

        DiffModSat[iFrame,:,:] = TmpMatrix

    ActualMonth = req.month-1
    ChlDiffSquared[ActualMonth,:,:] = np.nanmean(DiffModSat,0)
    mnan = np.isnan(ChlDiffSquared[ActualMonth,:,:])
    ChlDiffSquared[ActualMonth,mnan] = fillValue
#for month in range(0,12):
    fname = OUTDIR + "/varErr.%02d.nc"%(req.month)
      
    Sat.dumpGenericNativefile(fname,ChlDiffSquared[ActualMonth,:,:],varname='variance',mesh=MyMesh)

print "\n\tDone :)\n"


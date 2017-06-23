import numpy as np
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.timerequestors import Clim_month
import Sat.SatManager as Sat
import netCDF4

def VarianceDiffModelObs(MyMesh, ModDir, SatDir, OutDir):
  """
  Computes the monthly variance of the difference between Model output and Satellite observations
  
  MyMesh must be consistent with both the Model data inputs (within ModDir)
  and the Satellite observations  (within SatDir)

  """

  # TimeInterval MUST BE coherent with the inputs folder
  # In particular, it have to coincide with the "smallest" time interval
  # among the two folders (ModDir and SatDir)
  Timestart="19990103"
  Time__end="20150828"

  TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
  TLModel = TimeList.fromfilenames(TI, ModDir,"*.nc",prefix='chl.',dateformat="%Y%m%d")
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

  for req in MONTHLY_reqs:
    iiMod, wMod = TLModel.select(req)
    iiSat, wSat = TLSat.select(req)

    assert(len(iiMod) == len(iiSat)), "TimeLists not coherent"

    nFiles = len(iiMod)
    DiffModSat = np.zeros((nFiles,jpj,jpi))

    for iFrame, j in enumerate(iiMod):

      # decomment the following line to make the script working
      # if j==0:
      #   continue

      ModFile = TLModel.filelist[j]
      SatFile = TLSat.filelist[j]

      ModChl = Sat.readfromfile(ModFile,'lchlm')
      SatChl = Sat.readfromfile(SatFile,'CHL')

      MyMask = np.where(SatChl == fillValue)

      # debug output
      # print "TLModel.Timelist[j] ", TLModel.Timelist[j], "  TLSat.Timelist[j] ",TLSat.Timelist[j], " j = ", j

      if j > 0:
        assert(TLModel.Timelist[j]==TLSat.Timelist[j])
      assert(DiffModSat[0,:,:].shape == ModChl.shape)
      assert(DiffModSat[0,:,:].shape == SatChl.shape)

      TmpMatrix = ModChl - SatChl
      TmpMatrix = TmpMatrix*TmpMatrix      
      TmpMatrix[MyMask] = fillValue

      DiffModSat[iFrame,:,:] = TmpMatrix

    ActualMonth = req.month-1
    ChlDiffSquared[ActualMonth,:,:] = Sat.averager(DiffModSat)

  for month in range(0,12):
    if(OutDir[-1] == '/'):
      fname = "varErr.%02d.nc"%(month+1)
    else:
      fname = "/varErr.%02d.nc"%(month+1)
      
    filename = OutDir+fname
    Sat.dumpGenericNativefile(filename,ChlDiffSquared[month,:,:],varname='variance',mesh=MyMesh)

  print "\n\tDone :)\n"

if __name__ == "__main__":
  
    from postproc.masks import Mesh24

    ModDir = "/pico/scratch/userexternal/ateruzzi/ChlSup_RACOAST_24/CHL_SUP24/"
    SatDir = "/pico/scratch/userexternal/pdicerbo/WorkDir/AveSat24/WeeklySat24Friday/"
    OUTDIR = "/pico/scratch/userexternal/pdicerbo/WorkDir/AveSat24/VarDiffModObs/"

    VarianceDiffModelObs(Mesh24, ModDir, SatDir, OUTDIR)

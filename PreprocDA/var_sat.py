import numpy as np
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.timerequestors import Clim_month
import Sat.SatManager as Sat
import netCDF4


def var_sat_CCI_10gg(limstd, dayF, MyMesh, INDIR, OUTDIR):
  """
  Computes the Monthly surface variance of the 10-days averages
  The inputs file are within INDIR directory and the outputs will be placed
  within OUTDIR

  The mesh must be consistent with the data contained in INDIR
  """

  Timestart="19500101"
  Time__end="20500101"

  TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
  TLCheck = TimeList.fromfilenames(TI, INDIR,"*.nc",prefix='',dateformat='%Y%m%d')
  MONTHLY_reqs=[Clim_month(i) for i in range(1,13)]

  print "number of requests: ", len(MONTHLY_reqs), " inputFreq: ", TLCheck.inputFrequency

  jpi = MyMesh.jpi
  jpj = MyMesh.jpj

  Chl = np.zeros((12,jpj,jpi), dtype=float)
  ChlSquare = np.zeros((12,jpj,jpi), dtype=float)

  # get fillValue assuming that all files in the
  # folder have the same fillValue
  ncIN = netCDF4.Dataset(TLCheck.filelist[0],'r')
  varObj = ncIN.variables['CHL']
  fillValue = varObj.getncattr('_FillValue')
  ncIN.close()

  for req in MONTHLY_reqs:
    ii, w = TLCheck.select(req)
  
    nFiles = len(ii)
    M  = np.zeros((nFiles,jpj,jpi),np.float32)
    M2 = np.zeros((nFiles,jpj,jpi),np.float32)
    
    for iFrame, j in enumerate(ii):
        inputfile = TLCheck.filelist[j]
        CHL = Sat.readfromfile(inputfile,'CHL')
        TmpMat = np.zeros(CHL.shape)
        M[iFrame,:,:]  = CHL

        MyMask = np.where(CHL==fillValue)
        TmpMat = CHL*CHL
        TmpMat[MyMask] = fillValue
        M2[iFrame,:,:] = TmpMat
    
    MonthIndex = req.month-1
    Chl[MonthIndex,:,:] = Sat.WeightedAverager(M,w)
    ChlSquare[MonthIndex,:,:] = Sat.WeightedAverager(M2,w)

  for ii in range(0,12):
    # computing variance
    Var2D = ChlSquare[ii,:,:] - Chl[ii,:,:]*Chl[ii,:,:]
    
    # filling missing points
    MyMask = np.where(ChlSquare[ii,:,:] == fillValue)
    Var2D[MyMask]  = fillValue
    
    # saving results
    MonthStr = '%02d'%(ii+1)
    fname = 'var2Dsat.CCI.F'+str(dayF)+'.'+str(limstd)+'.'+MonthStr+".nc"
    filename = OUTDIR+fname
    print "\tsaving ", fname
    Sat.dumpGenericNativefile(filename,Var2D,varname='variance',mesh=MyMesh)

  print "\n\tDone :)\n"

if __name__ == "__main__":

    from postproc.masks import Mesh24

    INDIR  = "/pico/scratch/userexternal/pdicerbo/WorkDir/AveSat24/Checked_10Days_SatInterp24/"
    OUTDIR = "/pico/scratch/userexternal/pdicerbo/WorkDir/AveSat24/VarSat10Days/"
    limstd = 2
    dayF   = 10

    var_sat_CCI_10gg(limstd, dayF, Mesh24, INDIR, OUTDIR)

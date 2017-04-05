import numpy as np
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval

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
  MONTHLY_reqs=TLCheck.getMonthlist()

  jpi = MyMesh.jpi
  jpj = MyMesh.jpj

  for req in MONTHLY_reqs:  
    print req, "\tYear:", req.year, "\tMonth", req.month 

  print "\n\tDone :)\n"

if __name__ == "__main__":

    from postproc.masks import Mesh24
    # print "pippo"

    INDIR  = "/pico/scratch/userexternal/pdicerbo/WorkDir/AveSat24/Checked_10Days_SatInterp24"
    OUTDIR = "/pico/scratch/userexternal/pdicerbo/WorkDir/AveSat24/VarSat10Days/"
    limstd = 2
    dayF   = 10

    var_sat_CCI_10gg(limstd, dayF, Mesh24, INDIR, OUTDIR)

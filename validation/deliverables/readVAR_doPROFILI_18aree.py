#  readVAR_doPROFILI18aree.py 
# read netcdf file and produce profile for selected 18 aree of 4 x4 degree
# one profile for each ave file.
# 
# the following path have to be set:
#
# export OPA_HOME=S:qO-03
# export    MASKFILE=$CINECA_SCRATCH/$OPA_HOME/wrkdir/MODEL/meshmask.nc
# export SUBMASKFILE=$CINECA_SCRATCH/$OPA_HOME/wrkdir/MODEL/submask.nc


import scipy.io.netcdf as NC
import glob
import os,sys
import numpy as np
curdir='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER'
# Definition of 18 areas: lower left point 
#LAT1=np.array[  34 35 36 40 36.0 41 34.0 39 42 32.5 36.5 40.5 33 31.0 35 39 32.0 32]
LAT1=np.array([34., 35., 36., 40., 36.0, 41., 34.0, 39., 42., 32.5, 36.5, 40.5, 33, 31.0, 35., 39., 32.0, 32.])
LON1=np.array([-5.5, -1., 03., 03., 07.5, 06., 11.5, 10., 12., 16.0, 17., 16.0, 20., 24.0, 24., 23., 28.0, 32.,])

#sys.path.append("/pico/scratch/userexternal/gbolzon0/Carbonatic/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP")
#os.chdir('/pico/scratch/userexternal/gbolzon0/Carbonatic/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP')

# read MESHMASK
maskfile='/pico/scratch/userexternal/gbolzon0/Carbonatic-02/wrkdir/MODEL/meshmask.nc'
M=NC.netcdf_file(maskfile,"r")
jpi=M.dimensions['x']
jpj=M.dimensions['y']
jpk=M.dimensions['z']
tmask   = (M.variables['tmask'].data[0,:,:,:]).astype(np.bool)
nav_lev =  M.variables['nav_lev'].data
Lon     =  M.variables['glamt'].data[0,0,:,:]
Lat     =  M.variables['gphit'].data[0,0,:,:]
area    = (M.variables['e1t'].data[0,0,:,:])*(M.variables['e2t'].data[0,0,:,:])
e3t     =  M.variables['e3t'].data[0,:,0,0]
M.close()
# build mask 1 false for area with depth> 200
mask200=tmask[25,:,:].copy()
# mask200 x number = number (mask200=True), = 0 (mask200=False) 
mask0_200=tmask[0,:,:]*~mask200

# read the monthly files and build the 6 layer average and performe the 1x1 degree average
#INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-08/wrkdir/MODEL/AVE_FREQ_2/'
INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-10/wrkdir/MODEL/AVE_FREQ_2/'
INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-11/wrkdir/MODEL/AVE_FREQ_2/'
INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-12/wrkdir/MODEL/AVE_FREQ_2/'
INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-17/wrkdir/MODEL/AVE_FREQ_2/'
#INPUT_PHYSDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-11/wrkdir/MODEL/AVE_PHYS/'
#INPUT_PHYSDIR='/pico/scratch/userexternal/gcossari/COPERNICUS/tmp_avephys_run_Carbonatic/' 

OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_11/'
OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_12/'
OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_17/'

VARlist = ['AC_']
#VARlist = ['PH_']
VARlist = ['pCO']
#VARlist = ['DIC']
#VARlist = ['F05'] # NET GROWTH OF P2
#VARlist = ['F11'] # fDIC due to calcifier
#VARlist = ['F12'] # fALK due to calcifier
VARlist = ['O3h'] 
VARlist = ['O3c'] 

#VARlist = ['ppn']  # mgC/m3/d

#VARlist = ['N1p']
#VARlist = ['N3n']
#
#VARlist = ['votemper']
#VARlist = ['vosaline']

#VARlist = ['P1c']
#VARlist = ['P2c']
#VARlist = ['P3c']
#VARlist = ['P4c']

#VARlist = ['N5s']
#VARlist = ['R6c']
#VARlist = ['R2c']

#VARlist = ['Z3c']
#VARlist = ['Z4c']
#VARlist = ['Z5c']
#VARlist = ['Z6c']
#VARlist = ['O2o']
#VARlist = ['B1c']

#VARlist = ['P2n']
#VARlist = ['P2p']
#VARlist = ['N4n']

#VARlist = ['F13'] # contain EmPmR for DIC and ALK


for varname in VARlist:
   print ("start reading " + varname )
   if varname == 'votemper' or  varname == 'vosaline':
      aveLIST = glob.glob(INPUT_PHYSDIR + "ave.*.phys.nc")
   else:
      aveLIST = glob.glob(INPUT_AVEDIR+"ave*."+ varname + ".nc")
   aveLIST.sort()
   it=0
   Lnull=10
   VARiz=np.zeros((15,jpk-Lnull,18),dtype=np.double)
 #  aveLIST[1:]=[] # DA USARE PER FAR FARE UN GIRO SOLO AL CICLO SUL TEMPO
   for avefile in aveLIST:   #[rank::nranks]:
      print("reading " + avefile + " ...")
      ncIN = NC.netcdf_file(avefile,"r")
      A = ncIN.variables[varname].data[0,:,:,:]
      VAR=A.copy()
      VAR[VAR>1e+19]=np.NaN # NaN for no values
      ncIN.close()
      for jk in range(0,jpk-Lnull):
        VAR1=VAR[jk,:,:]*mask200 
        VAR1[VAR1==0]=np.NaN
        for iz in range(0,18):
               a=VAR1[(Lon>=LON1[iz]) & (Lon<LON1[iz]+4) & (Lat>=LAT1[iz]) & (Lat<LAT1[iz]+4) ]
               VARiz[it,jk,iz]=np.nansum(a) / np.sum(np.isfinite(a))
      it=it+1
   VARiz[np.isnan(VARiz)]=1.E20
#  scritttura del file NC\
   nomefile = 'PROF_18aree_' + varname +'.nc'
   outfile = OUTPUT_DIR + nomefile
   ncOUT=NC.netcdf_file(outfile,"w")
# write dimensions
   ncOUT.createDimension('depth',jpk-Lnull)
   ncOUT.createDimension('area',18)
   ncOUT.createDimension('time',15)
# write variables for dimensions
   time = ncOUT.createVariable('time','f',('time',))
   time[:]=[1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.]
   time.units = 'months starting from 04/2014'

   depth = ncOUT.createVariable('depth','f',('depth',))
   depth[:] = nav_lev[0:jpk-Lnull]
   depth.units = 'depth of centre of the cells:nav_lev' 

   area = ncOUT.createVariable('area','f',('area',))
   area[:]=[1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.]
   area.units = '18 areas, see readme for definition of area'
   ncOUT.description = 'chain run with DA from 04-14 to 05-15 monthly: 18 areas x vertical profile'

#   save variable 
   data = ncOUT.createVariable(varname,'f',('time','depth','area',))
   data[:]=VARiz
   setattr(data,"missing_value",1.0e+20)
   setattr(data,"actual_range",[np.min(VARiz),np.max(VARiz)])
   ncOUT.close()








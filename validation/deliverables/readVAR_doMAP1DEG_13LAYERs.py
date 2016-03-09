#  readVAR_doMAP1DEG_LAYERs.py 
# read netcdf file and produce monthly value for selected layer: ...
# number of layer 13: 50m from 0 to 200, 200-500 and then layer of 500
# and a degriding to 1 x 1 deg 
# thus this maps are to be compared to the CARBSYS climatology
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

LON=np.arange(-8.5,35.501,1.0);  # equivalente del matlab -8.5:1:35.5
Xlon=np.arange(-8,35.001,1.0);
LAT=np.arange(29.5,46.501,1.0);  # equivalente del matlab 29.5:1:46.5
Ylat=np.arange(30,46.001,1.0);

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

dZ     =np.ones ((jpk,jpj,jpi),dtype=np.double)

for i in range(jpk):
    dZ[i,:,:]     = e3t[i] * tmask[i,:,:]


# buold the height of the six layers
DEPTH1=np.array([00, 050, 100, 150, 200, 0500, 1000, 1500, 2000, 2500, 3000, 3500, 4000])
DEPTH2=np.array([50, 100, 150, 200, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500])

# see excel file PROF_CELL_MESHGRID_1_16.xlxs

LEV1=np.array([00, 11, 17, 22, 25, 37, 46, 52, 57, 60, 63, 65, 67]) # levels of BFM meshmask MED 16
LEV2=np.array([10, 16, 21, 24, 36, 45, 51, 56, 59, 62, 64, 66, 68])

VARlev =np.zeros ((jpj,jpi),dtype=np.double)
VARlev1x1=np.zeros((15,13,np.size(Ylat),np.size(Xlon)),dtype=np.double)

if 1==2:   # non lo rifaccio ad ogni exe
 Hlev   =np.zeros ((13,jpj,jpi),dtype=np.double)
 mask1x1=np.zeros((13,np.size(Ylat),np.size(Xlon)),dtype=np.double)
 for ilev in range(0,13):
    print(ilev)
    Hlev[ilev,:,:]=np.sum(dZ[LEV1[ilev]:LEV2[ilev]+1,:,:],0)
    Hlevmask=Hlev[ilev,:,:].copy()
    Hlevmask[Hlevmask>0]=1.;
    for ilon in range(np.size(LON)-1):
       for ilat in range(np.size(LAT)-1):
          a=Hlevmask[(Lon>=LON[ilon]) & (Lon<LON[ilon+1]) & (Lat>=LAT[ilat]) & (Lat<LAT[ilat+1]) ]
        #  if np.sum(a)>160.:  # 60% of 256
          if np.sum(a)>130.:  # 30% of 256
             mask1x1[ilev,ilat,ilon]=1.
 np.save('PYfile_Hlev_Hlevmask_mask1x1.npy', [Hlev,Hlevmask,mask1x1])
else:
 [Hlev,Hlevmask,mask1x1]=np.load('PYfile_Hlev_Hlevmask_mask1x1.npy')
#


# read the monthly files and build the 6 layer average and performe the 1x1 degree average
#INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP'
#INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-01/wrkdir/MODEL/AVE_FREQ_2/'
#INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-02/wrkdir/MODEL/AVE_FREQ_2/'
#INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-03/wrkdir/MODEL/AVE_FREQ_2/'
#INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-04/wrkdir/MODEL/AVE_FREQ_2/'
#INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-05/wrkdir/MODEL/AVE_FREQ_2/'
#INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-06/wrkdir/MODEL/AVE_FREQ_2/'
#INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-07/wrkdir/MODEL/AVE_FREQ_2/'
#INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-08/wrkdir/MODEL/AVE_FREQ_2/'
INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-09/wrkdir/MODEL/AVE_FREQ_2/'
INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-10/wrkdir/MODEL/AVE_FREQ_2/'
#INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-11/wrkdir/MODEL/AVE_FREQ_2/'
#INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-12/wrkdir/MODEL/AVE_FREQ_2/'
INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-14/wrkdir/MODEL/AVE_FREQ_2/'
INPUT_AVEDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic-17/wrkdir/MODEL/AVE_FREQ_2/'
#INPUT_PHYSDIR='/pico/scratch/userexternal/gbolzon0/Carbonatic/wrkdir/MODEL/AVE_PHYS/'
#INPUT_PHYSDIR='/pico/scratch/userexternal/gcossari/COPERNICUS/tmp_avephys_run_Carbonatic/' 

#OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_01/'
#OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_02/'
#OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_03/'
#OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_04/'
#OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_05/'
#OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_06/'
#OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_07/'
#OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_08/'
OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_09/'
OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_10/'
OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_14/'
#OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_11/'
#OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_12/'
OUTPUT_DIR='/pico/home/userexternal/gcossari/COPERNICUS/CALCIFIER/run_carbonatico_17/'

VARlist = ['AC_']
VARlist = ['PH_']
#VARlist = ['pCO']
#VARlist = ['DIC']
#VARlist = ['F05'] # NET GROWTH OF P2
#VARlist = ['F11'] # fDIC due to calcifier
#VARlist = ['F12'] # fALK due to calcifier
#VARlist = ['O3h'] 
#VARlist = ['O3c'] 

#VARlist = ['ppn']  # mgC/m3/d

VARlist = ['N1p']
VARlist = ['N3n']
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
 #  aveLIST[1:]=[] # DA USARE PER FAR FARE UN GIRO SOLO AL CICLO SUL TEMPO
   for avefile in aveLIST:   #[rank::nranks]:
      print("reading " + avefile + " ...")
      ncIN = NC.netcdf_file(avefile,"r")
      A = ncIN.variables[varname].data[0,:,:,:]
      VAR=A.copy()
      VAR[VAR>1e+19]=np.NaN # NaN for no values
      ncIN.close()
      
      VARdZ=VAR*dZ # multiply VAR for height of the cells
      for ilev in range(0,13):  # VARdZ contains NaN therefore nansum is used
          VARlev[:,:]=np.nansum(VARdZ[LEV1[ilev]:LEV2[ilev]+1,:,:],0) / Hlev[ilev,:,:]
          # averaged 1x 1 degree 
          for ilon in range(np.size(LON)-1):
           for ilat in range(np.size(LAT)-1):
             if mask1x1[ilev,ilat,ilon]>0:
               a=VARlev[(Lon>=LON[ilon]) & (Lon<LON[ilon+1]) & (Lat>=LAT[ilat]) & (Lat<LAT[ilat+1]) ]
               VARlev1x1[it,ilev,ilat,ilon]=np.nansum(a) / np.sum(np.isfinite(a))
             else:
               VARlev1x1[it,ilev,ilat,ilon]=np.nan          
      it=it+1
   VARlev1x1[np.isnan(VARlev1x1)]=1.E20
#  scritttura del file NC\
   nomefile = 'MAP1x1_13lev_15m_' + varname +'.nc'
   outfile = OUTPUT_DIR + nomefile
   ncOUT=NC.netcdf_file(outfile,"w")
# write dimensions
   ncOUT.createDimension('depth',13)
   ncOUT.createDimension('lat',17)
   ncOUT.createDimension('lon',44)
   ncOUT.createDimension('time',15)
# write variables for dimensions
   time = ncOUT.createVariable('time','f',('time',))
   time[:]=[1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.]
   time.units = 'months starting from 04/2014'

   depth = ncOUT.createVariable('depth','f',('depth',))
   depth[:] = [21.2,71.2,124.3, 177.2, 335.9, 721.5, 1210.7, 1758.9, 2293.9, 2798.4, 3294.5, 3752.3, 4269.3]
   depth.units = 'depth of centre of the layers: 50m from 0 to 200, 200-500, 500m from 500 to 4500' 

   lon = ncOUT.createVariable('lon','f',('lon',))
   lat = ncOUT.createVariable('lat','f',('lat',))
   lon[:]=Xlon
   lat[:]=Ylat
   ncOUT.description = 'chain run with DA from 04-14 to 05-15 monthly: 1x1 degree , 13 layers'

#   save variable 
   data = ncOUT.createVariable(varname,'f',('time','depth','lat','lon',))
   data[:]=VARlev1x1
   setattr(data,"missing_value",1.0e+20)
   setattr(data,"actual_range",[np.min(VARlev1x1),np.max(VARlev1x1)])
   ncOUT.close()








import sys
import argparse
from commons.utils import addsep
import numpy as np
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.mask import Mask
from commons.layer import Layer
from commons.utils import writetable
import netCDF4 as NC
import commons.netcdf4 as NC4 
import matplotlib.pyplot as pl
from layer_integral.mapplot import mapplot,pl
from layer_integral import coastline

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generator of png of comparisons between CLIM and MODEL degraded to 1x1 res
    Generator of a table with BIAS, RMSD and CORR for the two compared dataset
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--inputdir', '-i',
                              type = str,
                              default = "/galileo/home/userexternal/lfeudale/matlab_MAP_1x1/2019/",
                              required = True,
                              help = ''' Path of model 1x1 ''')

    parser.add_argument(   '--climdir', '-c',
                              type = str,
                              default = "/galileo/home/userexternal/lfeudale/matlab_MAP_1x1/MAP1x1/",
                              required = True,
                              help = ''' Path of clim 1x1 ''')

    parser.add_argument(   '--outdir', '-o',
                              type = str,
                              default = None,
                              required = True,
                              help = "")

    return parser.parse_args()

args = argument()



font_size=15

clon,clat = coastline.get()
Maskfile="/galileo/home/userexternal/lfeudale/matlab_MAP_1x1/mask_1x1_8lev.nc"
#cldir = "/galileo/home/userexternal/lfeudale/matlab_MAP_1x1/MAP1x1/"
#moddir= "/galileo/home/userexternal/lfeudale/matlab_MAP_1x1/2019/"
cldir = addsep(args.climdir)
moddir = addsep(args.inputdir)
outdir = addsep(args.outdir)

D=NC.Dataset(Maskfile,"r")  
mask=np.array(D.variables['tmask'])
 
mask=NC4.readfile(Maskfile,"tmask")
bmask=mask.astype(bool) 
kk=(mask==0)
nmask=mask.astype("float32")
nmask[kk]=np.nan

LAYERLIST=[Layer(0,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000), Layer(1000,2000)]
LAYERnames=["0-30","30-60","60-100","100-150","150-300","300-600","600-1000","1000-2000"]
VARLIST = ['P_l','N3n','O2o']
VARLIST=['N1p','N3n','O2o','ALK','DIC','pH','pCO2']
VARLIST=['N1p','N3n','O3c','O3h','DIC','Ac','PH','pCO','ppn','R6c','R2c'] #3c','O3h','Ac','PH','pCO','ppn','R6c','R2c'] #3c','O3h','Ac','PH','pCO','ppn','R6c','R2c'] #r32o','ALK','DIC','pH','pCO2']
ERRLIM=[[-0.1,0.1],[-0.1,0.1],[-1000,1000],[-80,80],[-80,80],[-80,80],[-0.05,0.05],[-0.05,0.05],[-5,5],[-0.1,0.1],[-0.1,0.1]]
VARLIST=['N1p','N3n','O3c','O3h','DIC','Ac','PH']
VARLIST_clim=["PHOS","NIT","O3C","O3H","DIC","ALK","PH_TOT_T_Press_ins"]
VARUNIT=['\mumol/m^4','\mumol/m^3','mgC/m^3','mmol/m^3','\mumol/kg','\mumol/kg',' ',' ','mgC/m3/d','mgC/m^3','mgC/m^3']
clim=[[0,0.25],[0,6],[26000,29000],[2400,2800],[2050, 2400],[2380,2750],[8.0,8.2],[300,420],[-1,10],[0,12],[0,30]]
nLat=17
nLon=44

Xlon=np.arange(-8,36)  # ATLANTIC BUFFER FROM -5
xs=3
Ylat=np.arange(30,47)

#Xlon=[-8:1:35];Ylat=[30:1:46];
VARLIST=['N1p','N3n','O3c','O3h','DIC','Ac','PH','pCO','ppn','R6c','R2c']
# THE FOLLOWING IS WHAT IS AVAILABLE:
VARLIST=['N1p','N3n','DIC','ALK','pH'] #,'pCO']
VARLIST_clim=["PHOS","NIT","DIC","Ac","PH"] # "PH_TOT_T_Press_ins"]
VARUNIT=['\mumol/m^4','\mumol/m^3','\mumol/kg','\mumol/kg',' ']
clim=[[0,0.25],[0,6],[2050, 2400],[2380,2750],[8.0,8.2]] 
ERRLIM=[[-0.1,0.1],[-0.1,0.1],[-80,80],[-80,80],[-0.05,0.05]]
nLayers=len(LAYERLIST)
nSeasons=5 # nomestag ={'YEAR','WIN','SPR','SUM','AUT'}; stagini=[01 01 04 07 10]; stagfin=[12 03 06 09 12];
nomestag ={'YEAR','WIN','SPR','SUM','AUT'};
nStat = 6 # stat: mean, std, #, median, p25, p75)
Z_clim = np.array([-(l.bottom+l.top)/3  for l in LAYERLIST])

CLI=np.zeros((nStat,nSeasons,nLayers,nLat,nLon),np.float32)*np.nan
MOD=np.zeros((12,nLayers,nLat,nLon),np.float32)*np.nan


#MOD[MOD>1.0e+15]==np.nan
#pippo=MOD[MOD<1.0e+20]


for ivar,VARname in enumerate(VARLIST):
    nomeMODfile=moddir + 'MAP1x1_13lev_'+  VARname + '.nc'
    if ( VARname == "pH") : 
        nomeCLIfile=cldir+  "PH_TOT_T_Press_ins" + '_lon_lat_1GRADO_8depth_5stag.nc'
    else:
        nomeCLIfile=cldir+  VARLIST_clim[ivar] + '_lon_lat_1GRADO_8depth_5stag.nc'
    print VARname
    print nomeCLIfile
    if (( VARname == "ALK") | (VARname == "pH")) : 
        CLI=NC4.readfile(nomeCLIfile,VARLIST_clim[ivar])
    else:
        CLI=NC4.readfile(nomeCLIfile,VARname)
    MOD=NC4.readfile(nomeMODfile,VARname)

    # Change the 1.e+20 in Nan
    ii=(CLI==1.e+20)
    CLI[ii]=np.nan

    ij=MOD==1.e+20
    BIAS=np.zeros((nLayers),np.float32)*np.nan
    RMSD=np.zeros((nLayers),np.float32)*np.nan
    CORR=np.zeros((nLayers),np.float32)*np.nan
    STAT=np.zeros(((nLayers),3),np.float32)*np.nan

    for zlev in range(nLayers):
        print LAYERLIST[zlev]
        MOD[ij]=np.nan
#        MOD_Y=np.nanmean(MOD[:,zlev,:,xs:],axis=0) # STARTING FROM -6 lon instead of -8
#        CLI_Y=CLI[0,0,zlev,:,xs:]                  # STARTING FROM -6 lon instead of -8
        MOD_Y=np.nanmean(MOD[:,zlev,:,xs:]*nmask[zlev,:,xs:],axis=0) # STARTING FROM -6 lon instead of -8
        CLI_Y=CLI[0,0,zlev,:,xs:]*nmask[zlev,:,xs:]                  # STARTING FROM -6 lon instead of -8
#        DIF=np.nanmean(MOD[:,zlev,:,:],axis=0)-CLI[0,0,zlev,:,:]
        DIF=MOD_Y-CLI_Y
        fig , ax = pl.subplots(3,1,figsize=(10,15),dpi=150,facecolor='w',edgecolor='k')
        ax = ax.ravel()
        im0=ax[0].pcolormesh(Xlon[xs:], Ylat, MOD_Y,vmin=clim[ivar][0],vmax=clim[ivar][1])
        cb0=fig.colorbar(im0,ax=ax[0])
        cb0.ax.set_yticklabels(cb0.ax.get_yticklabels(), fontsize=font_size)
        ax[0].plot(clon,clat,'k')
        ax[0].tick_params(labelbottom=False)
        tit1= VARname
        ax[0].set_title(VARname + " model ave Level " + LAYERnames[zlev])
  
        im1=ax[1].pcolormesh(Xlon[xs:], Ylat, CLI_Y,vmin=clim[ivar][0],vmax=clim[ivar][1])
        cb1=fig.colorbar(im1,ax=ax[1])
        cb1.ax.set_yticklabels(cb1.ax.get_yticklabels(), fontsize=font_size)
        ax[1].plot(clon,clat,'k')
        ax[1].set_title("OBS")
        ax[1].tick_params(labelbottom=False)
        font_size=15
        im2=ax[2].pcolormesh(Xlon[xs:], Ylat, DIF[:,:],vmin=ERRLIM[ivar][0],vmax=ERRLIM[ivar][1],cmap="bwr")
        cb2=fig.colorbar(im2,ax=ax[2])
        cb2.ax.set_yticklabels(cb2.ax.get_yticklabels(), fontsize=font_size)
        ax[2].tick_params(labelsize=font_size)
        ax[2].plot(clon,clat,'k')
        ax[2].set_title("MOD-OBS")
        for k in [0,1,2]:
            ax[k].tick_params(labelsize=font_size)
            ax[k].set_xlim([-6,36])
            ax[k].set_ylim([30,46])
#        fig.show()
        STAT[zlev,0]=np.nanmean(DIF)
        STAT[zlev,1]=np.sqrt(np.nanmean(DIF**2))
        
        A=np.resize(MOD_Y,((MOD_Y).shape[0],(MOD_Y).shape[1]))
        B=np.resize(CLI_Y,((CLI_Y).shape[0],(CLI_Y).shape[1]))
        CORR=np.corrcoef(A[~np.isnan(A*B)],B[~np.isnan(A*B)]) #CLIM has only few gridpoints availables
        STAT[zlev,2]=CORR[0,1]
        fig.savefig("CFR_M_O_MAP1x1_" + VARname + "_lev_" + LAYERnames[zlev] + "_aveYEAR.png")
   
    rows_names_list=LAYERnames
    column_names_list=["BIAS","RMSD","CORR"]

    outfile= VARname + "-LAYER-Y-CLASS4-CLIM-BIAS-RMSD-CORR.txt"
    writetable(outfile, STAT, rows_names_list, column_names_list) #,fmt=("%8.2f\t   %8.2f\t   %8.2f\t   %d"))

#    sys.exit()

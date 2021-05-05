import numpy as np
import os,sys
from commons import netcdf4
from commons.utils import addsep
from commons.mask import Mask
from mhelpers.linear_shift import linear_shift
import seawater as sw

woa_dir="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/WOA2018/Med/"
WOA_DIR=addsep(os.getenv("WOA_DIR", woa_dir))
if not os.path.exists(WOA_DIR):
    print WOA_DIR
    raise ValueError("Environment variable WOA_DIR must be defined")

maskfile = WOA_DIR + "meshmask.nc"

Mask_WOA = Mask(maskfile)
levwoa = Mask_WOA.zlevels
nwoa1000   = (levwoa<=1000).sum()
ind_woa600 = (levwoa<= 600).sum()


def read_climatology_nitrate():

    N3n = netcdf4.readfile(WOA_DIR + "ave.20000101-00:00:00.N3n.nc", "N3n")[0]
    temp = netcdf4.readfile(WOA_DIR + "ave.20000101-00:00:00.votemper.nc", "votemper")[0]
    sali = netcdf4.readfile(WOA_DIR + "ave.20000101-00:00:00.vosaline.nc", "vosaline")[0]
    good = N3n<1.e+30
    jpk, jpj, jpi = Mask_WOA.shape
    Pres3D =np.zeros((jpk,jpj,jpi),np.float32)
    density=np.zeros((jpk,jpj,jpi),np.float32)
    for i in range(jpi):
        for j in range(jpj):
            Pres3D[:,j,i] = levwoa

    density[good] = sw.dens(sali[good],temp[good],Pres3D[good])
    N3n = N3n * density/1000.
    N3n[~good] = np.nan

    return N3n

woa3D = read_climatology_nitrate()



def woa_nitrate_correction(p):
    '''
    Calculates the nitrate value of WOA and correct by the bottom value
    Argument:
    * p * float profile object

    Return:
    * Pres  * numpy 1d array, corrected pressure
    * Value *
    * Qc    *
    '''

#             N3n_woa_mean = woa3D['N3n']['mn'][:nwoa1000,indxswoa[0],indxswoa[1]]
#             N3n_woa_obja = woa3D['N3n']['an'][:nwoa1000,indxswoa[0],indxswoa[1]]


# Find and extract WOA nitrate profile close to float coordinates
    jpi, jpj = Mask_WOA.convert_lon_lat_wetpoint_indices(p.lon,p.lat)

    N3n_WOA_prof_600_1000 = woa3D[ind_woa600:nwoa1000+1,jpj,jpi]
    N3n_WOA_bottom = np.nanmean(N3n_WOA_prof_600_1000)

# Find the "bottom" value of float nitrate [600m - 1000m mean]
    Np, N, Nqc = p.read('NITRATE',read_adjusted=True)
    ii = (Np>=600) & (Np<=1000) 
    N_Float_bottom=np.mean(N[ii])

    shift = N_Float_bottom - N3n_WOA_bottom

    New_profile = linear_shift(N,Np,shift)
    Nqc[:] = 8
    return Np, New_profile, Nqc, shift


if __name__ == "__main__":
    from instruments import bio_float
    from commons.time_interval import TimeInterval
    from commons.Timelist import TimeList
    import matplotlib.pyplot as pl

    from basins import V2 as OGS

    from commons.layer import Layer
    from static.climatology import get_climatology

    from woa_N3n import woa_nitrate_correction

    DATESTART = "20190101"
    DATE__END = "20191231"

    T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d')
    var="NITRATE"

    Profilelist=bio_float.FloatSelector(var,T_INT,OGS.med)
    print len(Profilelist)
    p=Profilelist[310]

    Pres, Prof, Qc= p.read(var,True)
    print len(Prof)

    from commons.layer import Layer
    from static.climatology import get_climatology

    from woa_N3n import woa_nitrate_correction

    DATESTART = "20190101"
    DATE__END = "20191231"

    if ((len(Prof)>5) & (Pres[-1]>600)):
    # PERFORM THE CORRECTION:
       Np, N, Nqc = woa_nitrate_correction(p)

       fig , ax = pl.subplots()
#    ax.plot(Prof,Pres,'b', N, Np,'r')
       ax.plot(Prof,Pres,'b',label="Original")
       ax.plot(N, Np,'r',label="Corrected [WOA]")


# Check with CLIMATOLOGY:

       SUBLIST = OGS.Pred.basin_list
       PresDOWN=np.array([0,25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500])
       LayerList=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]
       N3n_clim, N3n_std = get_climatology('N3n', SUBLIST, LayerList)
       z_clim = np.array([-(l.bottom+l.top)/2  for l in LayerList])

       for iSub,sub in enumerate(SUBLIST):
        if sub.is_inside(p.lon,p.lat):
           print N3n_clim[iSub,:]
           print sub.name
           sub_name=sub.name
           ax.plot(N3n_clim[iSub,:],-1*z_clim,'g',label="Clim")

       ax.invert_yaxis()
       ax.set_ylabel("depth $[m]$",color = 'k', fontsize=10)
       fig_title="Float ID " + np.str(p._my_float.wmo) + " NITRATE - " + p._my_float.time.strftime("%Y%m%d") + " [" + sub_name + "]"
       ax.set_title(fig_title)
       ax.legend()
#    fig.show()
       fig_name='FLOAT_' + np.str(p._my_float.wmo) + '_' + p._my_float.time.strftime("%Y%m%d") + '_' + sub_name  + '.png'
       fig.savefig(fig_name)

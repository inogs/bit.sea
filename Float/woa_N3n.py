import numpy as np
import os,sys
import netCDF4
from instruments import superfloat as bio_float
from commons.mask import Mask

MASKFILE_WOA="/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/GLO-ATL/WOA/meshmask.nc"

#DIRWOA=addsep(args.inwoa)
DIRWOA="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/WOA2018/Med/"
filewoa="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/WOA2018/Med/ave.20000101-00:00:00.N3n.nc"

Mask_WOA = Mask(MASKFILE_WOA)

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


# Reading WOA
    print '- Reading WOA -'

#woa3D = {}
#for var in WOAvar:
#    filewoa = DIRWOA + 
    print ' ... ' + filewoa
    dataset = netCDF4.Dataset(filewoa)

    lonwoa = dataset.variables['lon'][:].data
    latwoa = dataset.variables['lat'][:].data
    levwoa = dataset.variables['depth'][:].data
    levwoa1000 = levwoa[levwoa<=1000]
    nwoa1000 = len(levwoa1000)
    ind_woa600 = len(levwoa[levwoa<=600])
    firstreading = False

    woa3D = dataset.variables["n_mn"][0].data
    masknan = woa3D>10.e+30
    woa3D[masknan] = np.nan
    dataset.close()

#             N3n_woa_mean = woa3D['N3n']['mn'][:nwoa1000,indxswoa[0],indxswoa[1]]
#             N3n_woa_obja = woa3D['N3n']['an'][:nwoa1000,indxswoa[0],indxswoa[1]]


# Find and extract WOA nitrate profile close to float coordinates
    jpi, jpj = Mask_WOA.convert_lon_lat_wetpoint_indices(p.lon,p.lat)

    N3n_WOA_prof_600_1000 = woa3D[ind_woa600:nwoa1000+1,jpj,jpi]
    N3n_WOA_bottom = np.nanmean(N3n_WOA_prof_600_1000)

# Find the "bottom" value of float nitrate [600m - 1000m mean]
    Np, N, Nqc = p.read('NITRATE',True)
    ii = (Np>=600) & (Np<=1000) 
    N_Float_bottom=np.mean(N[ii])

    shift = N_Float_bottom - N3n_WOA_bottom

    New_N=N-shift

    iisurf = (New_N < 0) & (Np < 200)
    New_N[iisurf] = 0.05

# All the value profiles should be not negative
    ii = New_N > 0
    Np = Np[ii]
    New_N = New_N[ii]
    Nqc   =   Nqc[ii]

    return Np, New_N, Nqc


if __name__ == "__main__":

    from commons.time_interval import TimeInterval
    from commons.Timelist import TimeList
    import matplotlib.pyplot as pl
    import matplotlib.dates as mdates
    import numpy.ma as ma
    from basins import V2 as OGS

    from commons.layer import Layer
    from static.climatology import get_climatology

    from woa_N3n import woa_nitrate_correction

    DATESTART = "20190101"
    DATE__END = "20191231"

    T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d')
    var="NITRATE"

    Profilelist=bio_float.FloatSelector(var,T_INT,OGS.med)
    p=Profilelist[100]

    Pres, Prof, Qc= p.read(var,True)

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
           ax.plot(N3n_clim[iSub,:],-1*z_clim,'g',label="Clim")

    ax.invert_yaxis()
    ax.set_ylabel("depth $[m]$",color = 'k', fontsize=10)
    fig_title="Float ID " + np.str(p._my_float.wmo) + " NITRATE - " + p._my_float.time.strftime("%Y%m%d")
    ax.set_title(fig_title)
    ax.legend()
    fig.show()


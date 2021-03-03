import numpy as np
import os,sys
import netCDF4
from commons.utils import addsep
from commons.mask import Mask

woa_dir="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/WOA2018/Med/"
WOA_DIR=addsep(os.getenv("WOA_DIR", woa_dir))
if not os.path.exists(WOA_DIR):
    print WOA_DIR
    raise ValueError("Environment variable WOA_DIR must be defined")

maskfile = WOA_DIR + "meshmask.nc"
filewoa=   WOA_DIR + "ave.20000101-00:00:00.N3n.nc"

Mask_WOA = Mask(maskfile)



def read_climatology_nitrate():
    dataset = netCDF4.Dataset(filewoa)

    lonwoa = dataset.variables['lon'][:].data
    latwoa = dataset.variables['lat'][:].data
    levwoa = dataset.variables['depth'][:].data
    levwoa1000 = levwoa[levwoa<=1000]
    nwoa1000 = len(levwoa1000)
    ind_woa600 = len(levwoa[levwoa<=600])

    woa3D = dataset.variables["n_mn"][0].data
    masknan = woa3D>10.e+30
    woa3D[masknan] = np.nan
    dataset.close()
    return woa3D, ind_woa600, nwoa1000

woa3D, ind_woa600, nwoa1000 = read_climatology_nitrate()

def linear_shift(N_old,pres,qc,shift,p_bot=600):
    '''
    Perform a linar shift of nitrate profile along the vertical direction starting "for default" from 600m depth upward
    Argument:
    * N_old * Original nitrate profile
    * qc    * original QC
    * shift * shift calculated with WOA at bottom (mean 600-1000m)
    * p_bot * depth from which starts the shift (for default is fixed to 600m)

    Return:
    * New_profile * Nitrate profile with a linear shift applied
    * qc          * qc modified to 8 (that means interpolated value)
    '''

    N_new=N_old-shift

    P600=pres[(pres>=p_bot)][0]
    New_profile = N_new.copy()
    Shift_Surf = N_old[0]-max(0.05,N_new[0])
    for iz, zz in enumerate(pres[(pres<=p_bot)]):
        New_profile[iz] = N_old[iz] - (Shift_Surf + (shift - Shift_Surf)*(pres[iz]-pres[0])/(P600-pres[0]))
        New_profile[iz]=max(0.05,New_profile[iz])  # Eliminate possible negative values

    qc[:]=8

    return New_profile, qc


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

    New_profile, Nqc = linear_shift(N,Np,Nqc,shift)

    return Np, New_profile, Nqc


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

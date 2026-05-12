from typing import Tuple

import numpy as np
import gsw


def find_DCM(Chl_profile, zlev):
    A = Chl_profile
    CHL_surf = Chl_profile[0]
    A_filtered = A[A > 0.1]
    D_filtered = zlev[A > 0.1]
    A_fil_rev = A_filtered[::-1]
    D_fil_rev = D_filtered[::-1]

    if len(A_fil_rev) == 0:
        return np.nan, np.nan

    CM = A_fil_rev[0]
    DCM = D_fil_rev[0]
    for ip, chl in enumerate(A_fil_rev):
        if chl > CM:
            CM = chl
            DCM = D_fil_rev[ip]

    if DCM < 40:
        DCM = np.nan
    if CM / CHL_surf < 1.5:
        DCM = np.nan

    return CM, DCM


def smooth_profile(profile, N=5):
    profile_smoothed = convolve1d(profile, np.ones(N)/N, mode='reflect')
    return profile_smoothed

def MLD_sigma(Temperature, Salinity, Pres, insitu_T=True):
    '''
    Calculation of Mixed Layer Depth based on density difference
    of 0.03 kg m^-3 after Gali et al. 2022
    '''
    sigma_thr = 0.03 #[kg m^-3]
    N = 5 #running mean window
    nh = int(N/2)
    zref = 5.0 #[m]
    i5m = np.abs(Pres - zref).argmin()
    if insitu_T:
        # this not used cause gsw.pot_rho_t_exact() wants insitu_T
        #Temperature = ptmp(Salinity, Temperature, Pres)
        CT = gsw.CT_from_t(Salinity, Temperature, Pres)
    #
    # Density = pden(Salinity, Temperature, Pres)
    SA = gsw.SA_from_SP(Salinity, Pres, 0, 0)
    Density = gsw.pot_rho_t_exact(SA, Temperature, Pres, 10.1325) 
    d5m = np.mean(Density[i5m-nh:i5m+nh+1])
    for ip, p in enumerate(Density):
        abs_diff = p - d5m
        if abs_diff > sigma_thr:
            break
    MixedLayerDepth = np.min([1000.0, Pres[ip]])
    return MixedLayerDepth

def MLD(Temperature, Salinity, Pres, insitu_T=True, smooth=False):
    """
    Calculation of Mixed Layer Depth based on temperature difference of 0.2
    inputs are Temperature ([C], in-situ, optionally conservative temperature), Salinity ([psu]), Pressure ([kg/m^3])
    mld is defined as the level where the difference of temperature with respect the reference level of 10m is of 0.2C
    It resurns also DENSITY and POTENTIAL DENSITY at mld
    """
    th = 10  # threshold of depth minimum
    i_10 = np.abs(Pres - th).argmin()
    if smooth:
        # smooth T, S profiles
        N = 5
        Temperature = smooth_profile(Temperature, N)
        Salinity = smooth_profile(Salinity)
    SA = gsw.SA_from_SP(Salinity, Pres, 0, 0)
    if insitu_T:
        CT = gsw.CT_from_t(Salinity, Temperature, Pres)
    else:
        CT = Temperature
    # CONSIDER VALUES the reference level at 10m (de Boyer-Montegut 2004)
    D1000 = Pres[(Pres < 1000) & (Pres > th)]
    CT1000 = CT[(Pres < 1000) & (Pres > th)]
    SA1000 = SA[(Pres < 1000) & (Pres > th)]
    # in situ density anomaly [kg/m^3]
    Dens1000 = gsw.rho(SA1000, CT1000, D1000) - 1000  # in situ density anomaly [kg/m^3] DENSITY # SIGMA
    # potential density anomaly with reference pressure of 0 dbar
    PDens1000 = (gsw.sigma0(SA1000, CT1000))
    for ip, p in enumerate(CT1000):
        abs_diff = abs(p - CT[i_10])
        if abs_diff > 0.2:
            break
    MixedLayerDepth = D1000[ip]
    d_atMLD = Dens1000[ip]
    pd_atMLD = PDens1000[ip]
    return MixedLayerDepth, d_atMLD, pd_atMLD


def t_p_cline(
    Profile, Pres
):  # calculation of thermocline (Temp) - pycnocline (Dens)
    T = Profile
    dT = abs(np.diff(T)) / np.diff(Pres)
    ip = dT.argmax()
    return Pres[ip]


def StratIndex(
    BVF, TheMask, zdepth=1000
):  # BruntVaisalaFrequency is a 1D array
    iz1000 = TheMask.getDepthIndex(zdepth) + 1
    max_zindex = len(BVF)
    izmax = min(max_zindex, iz1000)
    SI = np.nan
    #        SI = np.nansum(BVF  * TheMask.zlevels[:izmax] * TheMask.dz[:izmax])/TheMask.dz[:izmax].sum()
    SI = np.nansum(
        BVF[:izmax, 0] * TheMask.zlevels[:izmax] * TheMask.dz[:izmax]
    )
    return SI


def find_WBL(Profile, Pres):  # Winter Bloom Layer (WBL)
    WBL = np.nan
    A = Profile
    A_filtered = A[Pres < 200]
    D_filtered = Pres[Pres < 200]
    for ip, p in enumerate(A_filtered):
        if p <= A[0] * 0.1:
            WBL = D_filtered[ip]
            break
    return WBL


def find_NITRICL(Profile, Pres):
    for ip, p in enumerate(Profile):
        if p >= 2:
            return Pres[ip]


def find_NITRICL_dz(Profile, Pres):
    dN = np.diff(Profile) / np.diff(Pres)
    for ip, p in enumerate(dN):
        if Pres[ip] > 40:
            if p >= 0.1:
                return Pres[ip]


def find_NITRICL_dz_max(Profile, Pres) -> Tuple[float, float]:
    """This is the Nitrcl2 used for the calcuation in the QUID.
    It can be used for nitracline and also pycnocline. Include
    also the variable value at that depth."""

    dN = np.diff(Profile) / np.diff(Pres)
    ip = dN.argmax()
    return Pres[ip], Profile[ip]


def find_OMZ(Profile, Pres) -> float:
    ii = (Pres > 200) & (Pres <= 1000)
    Pred = Profile[ii]
    j_minO2 = np.argmin(Pred)
    OMZ = Pres[j_minO2]
    return OMZ


def find_maxO2(Profile, Pres) -> float:
    Pred = Profile[Pres <= 200]
    j_maxO2 = np.argmax(Pred)
    MaxO2 = Pres[j_maxO2]
    return MaxO2


def find_bot_Nit(Profile, Pres) -> float:
    # Find the bottom value of Nitrate
    P600_800 = Profile[(Pres > 600) & (Pres <= 800)]
    Nit_bot = np.nanmean(P600_800)
    return Nit_bot

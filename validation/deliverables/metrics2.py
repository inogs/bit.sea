import numpy as np
import seawater as sw

def find_DCM(Chl_profile,zlev):

        A = Chl_profile
	CHL_surf = Chl_profile[0]
        A_filtered=A[A>0.1]
        D_filtered=zlev[A>0.1]
        A_fil_rev = A_filtered[::-1]
        D_fil_rev = D_filtered[::-1]

	if (len(A_fil_rev) == 0): 
	    return np.nan , np.nan

	CM  = A_fil_rev[0]
	DCM = D_fil_rev[0]
        for ip, chl in enumerate(A_fil_rev):
            if (chl > CM): 
                CM  = chl
                DCM = D_fil_rev[ip]

	if (DCM < 40): DCM = np.nan
	if (CM/CHL_surf < 1.5): DCM = np.nan

        return CM, DCM


def MLD(Temperature,Pres):
        ''' Calculation of Mixed Layer Depth based on temperature difference of 0.2
        mld is defined as 
        It resurns also DENSITY (SIGMA) and POTENTIAL DENSITY (SIGMA THETA)
        '''
        th=10 #threshold of depth minimum
        MLD = np.nan
        T = Temperature
        S = Salinity
        D1000=Pres[(Pres<1000) & (Pres>th)] # CONSIDER VALUES above "th"5m
        T1000=T[(Pres<1000) & (Pres>th)]
        S1000=S[(Pres<1000) & (Pres>th)]
        Dens1000=sw.dens(S1000,T1000,D1000)-1000 # DENSITY # SIGMA
        PDens1000=sw.pden(S1000,T1000,D1000)-1000 # POTENTIAL DENSITY # SIGMA THETA
        for ip,p in enumerate(T1000):
#            abs_diff=abs(p-T1000[0])
#            if abs_diff > 0.2:
            abs_diff=abs(p-np.mean(T[Pres<th]))
            if abs_diff > 0.2:
                        break
        MixedLayerDepth=D1000[ip]
        d_atMLD=Dens1000[ip]
        pd_atMLD=PDens1000[ip]
        return MixedLayerDepth, d_atMLD, pd_atMLD

def t_p_cline(Profile,Pres):  # calculation of thermocline (Temp) - pycnocline (Dens)
        T = Profile
        dT = np.diff(T)/np.diff(Pres)
        ip = dT.argmax()
        return Pres[ip] 

def StratIndex(BVF,TheMask): # BruntVaisalaFrequency is a 1D array
        iz1000 = TheMask.getDepthIndex(1000)+1
        max_zindex=len(BVF)
        izmax = min(max_zindex,iz1000)
        SI = np.nan
#        SI = np.nansum(BVF  * TheMask.zlevels[:izmax] * TheMask.dz[:izmax])/TheMask.dz[:izmax].sum()
        SI = np.nansum(BVF[:izmax,0]   * TheMask.zlevels[:izmax] * TheMask.dz[:izmax])
        return SI

def find_WBL(Profile,Pres):  # Winter Bloom Layer (WBL)
        WLB = np.nan
        A = Profile
	A_filtered=A[Pres<200]
	D_filtered=Pres[Pres<200]
        for ip, p in enumerate(A_filtered):
            if (p <= A[0]*0.1): 
                WLB = D_filtered[ip]
                Chl_min = A_filtered[ip]
	        break
        return WLB

def find_NITRICL(Profile,Pres):
        for ip, p in enumerate(Profile):
            if (p >= 2):
		return Pres[ip]

def find_NITRICL_dz(Profile,Pres):
         dN = np.diff(Profile)/np.diff(Pres)
         for ip, p in enumerate(dN):
          if ( Pres[ip] > 40 ):
            if (p >= 0.1):
                return Pres[ip]

def find_NITRICL_dz_max(Profile,Pres):
         ''' This is the Nitrcl2 used for the calcuation in the QUID.
             It can be used for nitracline and also pycnocline. Include 
             also the variable value at that depth. ''' 
 
         dN = np.diff(Profile)/np.diff(Pres)
         ip = dN.argmax()
         return Pres[ip], Profile[ip]    
 
def find_OMZ(Profile,Pres):
         ii = (Pres>200) & (Pres<=1000)
         Pred=Profile[ii]
         j_minO2=np.argmin(Pred)
         OMZ=Pres[j_minO2] 
         print OMZ
         return OMZ

def find_maxO2(Profile,Pres):
         Pred=Profile[Pres<=200]
         j_maxO2=np.argmax(Pred)
         MaxO2=Pres[j_maxO2]
         return MaxO2

def find_bot_Nit(Profile,Pres):
# Find the bottom value of nitrato 
         P600_800=Profile[(Pres>600) & (Pres<=800)]
         Nit_bot=np.nanmean(P600_800)
         return Nit_bot

import numpy as np

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


def find_DCM1(Chl_profile,zlev):
	CM  = np.nan
        DCM = np.nan
        ## type = 0 float or type = 1 model 
        A = Chl_profile
        A_filtered=A[A>0.1]
#       if (type == 0):
#               D_filtered=TheMask.zlevels[A[0,:]>0.1]
#       if (type == 1):
#               D_filtered=TheMask.zlevels[A[0,:]>0.1]
#       D_filtered=zlev[A[0,:]>0.1]
        D_filtered=zlev[A>0.1]
        A_fil_rev = A_filtered[::-1]
        D_fil_rev = D_filtered[::-1]

        for ip, p in enumerate(A_fil_rev[1:]):
            if (p < A_fil_rev[ip]) & (ip >0):
                CM  = A_fil_rev[ip]
                DCM = D_fil_rev[ip]
                break
        return (CM, DCM)

def find_DCM2(Chl_profile,zlev):
        CM  = np.nan
        DCM = np.nan

        A = Chl_profile
        A_filtered=A[A>0.1]
        D_filtered=zlev[A>0.1]
        A_fil_rev = A_filtered[::-1]
        D_fil_rev = D_filtered[::-1]
	d1 = np.diff(A_fil_rev,1)
        d2 = np.diff(A_fil_rev,2)
#        mindiff2=np.argmin(d2)
        d1sign = np.sign(d1) >= 0
	for ip, p in enumerate(d1):
	    if (ip > 0) and (d1sign[ip] == False) and (d2[ip-1] <0):
#	    if (d1[ip] == 0) and (d2[ip] <0):
#                CM  = A_fil_rev[ip]
#                DCM = D_fil_rev[ip]
                max_cand = (A_fil_rev[ip-1],A_fil_rev[ip],A_fil_rev[ip+1])
                d_max_cand = (D_fil_rev[ip-1],D_fil_rev[ip],D_fil_rev[ip+1])
		CM = max(max_cand)
		DCM = d_max_cand[max_cand.index(max(max_cand))] 
#		DCM = np.argmax(A_fil_rev[ip-1],A_fil_rev[ip],A_fil_rev[ip+1])
		break
        return (CM, DCM)

def find_MLD(Profile,Pres):
        MLD = np.nan
        A = Profile
#	A_filtered=A[Pres>20 and Pres<= 150]
#	D_filtered=Pres[Pres>20 and Pres<= 150]
	A_filtered=A[Pres<200]
	D_filtered=Pres[Pres<200]
	A_filtered=A[Pres<300]
	D_filtered=Pres[Pres<300]
#        A_filtered=A[0,TheMask.zlevels[:max_depth]>20]
#        D_filtered=TheMask.zlevels[TheMask.zlevels[:max_depth]>20]

#        for ip, p in enumerate(A_filtered[1:]):
#            if (p <= 0.05): 
        for ip, p in enumerate(A_filtered):
            if (p <= A[0]*0.1): 
# and D_filtered[ip] <= 150):
                MLD  = D_filtered[ip]
                Chl_min = A_filtered[ip]
	        break
#            if (MLD > 150):
#                MLD = np.nan
#                break
        return MLD

def find_NITRICL(Profile,Pres):
        nitrcl = np.nan
        ip = np.nan
        for ip, p in enumerate(Profile):
            if (p >= 2):
                nitrcl = Pres[ip]
                break
        return nitrcl,ip

def find_NITRICL_dz(Profile,Pres):
#        for ip, p in enumerate(Profile):
#         if ( Pres[ip] > 40 ):
#            dn = np.diff(p)/np.diff(Pres)
         dN = np.diff(Profile)/np.diff(Pres)
         for ip, p in enumerate(dN):
          if ( Pres[ip] > 40 ):
            if (p >= 0.1):
                return Pres[ip]

def find_NITRICL_dz_max(Profile,Pres):
         dN = np.diff(Profile)/np.diff(Pres)
         ip = dN.argmax()
         return Pres[ip]    

def find_NITRICL_dz_maxfondo(Profile,Pres):
        nitrcl = np.nan
        if Profile[-1] > 1.5:
            dN = np.diff(Profile)/np.diff(Pres)
            ip = dN.argmax()
            nitrcl = Pres[ip]
        return nitrcl

def find_THERMOCL(Profile,Pres,Temp):
        thermocl = np.nan
        ip = np.nan
        for ip, p in enumerate(Profile):
            if (p <= Temp):
                thermocl = Pres[ip]
                break
        return thermocl,ip

def find_SLOPE_dz_max(Profile,Pres):
         dN = np.diff(Profile)/np.diff(Pres)
         icl = dN.argmax()
         slope = np.nanmean(dN[icl-2:icl+2])
         return slope


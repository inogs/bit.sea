def find_DCM(Chl_profile,zlev):
	CM  = 0
        DCM = 0
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

def find_MLD(Profile,Pres):
        MLD = 0
        A = Profile
	A_filtered=A[Pres>20]
	D_filtered=Pres[Pres>20]
#        A_filtered=A[0,TheMask.zlevels[:max_depth]>20]
#        D_filtered=TheMask.zlevels[TheMask.zlevels[:max_depth]>20]

        for ip, p in enumerate(A_filtered[1:]):
            if (p <= 0.05): 
                MLD  = D_filtered[ip]
                Chl_min = A_filtered[ip]
                break
        return MLD

def find_NITRICL(Profile,Pres):
        for ip, p in enumerate(Profile):
            if (p >= 1):
		return Pres[ip]


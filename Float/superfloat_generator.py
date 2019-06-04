from instruments import bio_float
from instruments import lovbio_float
from commons.time_interval import TimeInterval
from basins.region import Rectangle
from commons import calculated_depths
import pylab as pl
import numpy as np
TI = TimeInterval('2015','2017','%Y')
R = Rectangle(-6,36,30,46)

def quenching(profile_obj, PresChl, Chl, chl_lov_zero):
    Quenched_Chl = Chl.copy()
    if profile_obj.time.month in range(1,13):
        PresT, Temp, _ = profile_obj.read('TEMP', read_adjusted=False)
        mld = calculated_depths.mld(Temp, PresT, zref=0)
        if mld > 200:mld=20

        ii=PresChl<=mld
        base_mld_level=ii.sum()
        nTries = 25
        if  nTries > len(PresChl)-base_mld_level-1:
            nTries = len(PresChl)-base_mld_level-1

        quench_value = np.zeros((nTries,), np.float32)
        diff_fromlov = np.zeros((nTries,), np.float32)
        for down_level in range(nTries):
            nLev = base_mld_level + down_level
            if nLev==0: nLev=1
            quench_value[down_level] = Chl[:nLev].max()
            diff_fromlov[down_level] = np.abs(quench_value[down_level]-chl_lov_zero)

        k = diff_fromlov.argmin()
        nLev= base_mld_level + k
        Quenched_Chl[:nLev] = quench_value[k]
    return Quenched_Chl

def are_identical(chl1, chl2):
    diff = chl1 - chl2
    if np.abs(diff).max()< 0.05:
        return True
    else:
        return False
def are_shifted(chl1, chl2):
    diff = chl1 - chl2
    l = np.median(diff) - np.mean(diff)
    if l<1.e-3:
        return True
    else:
        return False
def has_drift(Pres,Profile):
    if Pres.max()<200:
        return False
    i200 = np.abs(Pres-200).argmin()
    Delta = Profile[-1] - Profile[i200]
    if np.abs(Delta) > 0.05:
        return True
    else:
        return False


#PROFILES_COR =   bio_float.FloatSelector('CHLA', TI, R)
PROFILES_LOV =lovbio_float.FloatSelector('CHLA', TI, R)

for ip, pLov in enumerate(PROFILES_LOV[:]):
    pCor=bio_float.from_lov_profile(pLov)
    if pCor is not None:
        PresL, ValueL, QcL = pLov.read('CHLA',read_adjusted=True)
        if 'CHLA' in pCor.available_params:
            if pCor._my_float.status_var('CHLA')=='A':
                PresC, ValueC, QcC = pCor.read('CHLA',read_adjusted=True)
            else:
                print pCor._my_float.filename + " without CHLA_ADJUSTED"
                continue
        else:
            print pCor._my_float.filename + " without CHLA"
            continue

        if len(ValueC) < 5: continue
        if len(ValueL)<  1: continue
        first_derivative=np.diff(ValueC)/np.diff(PresC)
        max_jump = np.abs(np.diff(ValueC)).max()
        z_max = PresC[max_jump.argmax()]

        if (max_jump > 0.2) & (z_max > 100) :
            print "jump = ", max_jump, pCor._my_float.filename  + " will be replaced by LOV"
            continue
        if has_drift(PresC, ValueC):
            print pCor._my_float.filename + " has drift"
            continue
        CHL = quenching(pCor, PresC, ValueC,ValueL[0])

        ii=(PresC >= 400) & (PresC <= 600)
        if ii.sum() > 0:
            shift = CHL[ii].mean()
            CHL = CHL - shift

        if not are_identical(CHL, ValueL):
            if not are_shifted(CHL, ValueL):
                print ip
                fig,ax =pl.subplots()
                ax.plot(ValueL,PresL,'r', label="LOV")
                ax.plot(CHL  , PresC,'b.-', label="COR")
                ax.invert_yaxis()
                ax.grid()
                ax.legend()
                fig.show()
                pl.savefig(pCor.ID()+ ".png")
                #import sys
                #sys.exit()

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
    '''
    Arguments:
    * profile_obj  * a Coriolis profile object
    * Pres_chl     * numpy 1D-array of coriolis pressure
    * Chl          * numpy 1D-array of coriolis chl
    * chl_lov_zero * value of LOV quench; an algorithm will try to guess the Mixed Layer Depth
                     calculated by LOV, in order to use it here to calculate quench value
                     if None, the Coriols Mixed Layer Depth will be used.
    '''
    Quenched_Chl = Chl.copy()
    if profile_obj.time.month in range(1,13):
        PresT, Temp, _ = profile_obj.read('TEMP', read_adjusted=False)
        mld = calculated_depths.mld(Temp, PresT, zref=0)
        if mld > 200:mld=20

        ii=PresChl<=mld
        base_mld_level=ii.sum()
        if chl_lov_zero is not None:
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
        else:
            nLev = base_mld_level + 1
            quench_value = Chl[:nLev].max()
            Quenched_Chl[:nLev] = quench_value
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

def synthesis_profile(pLov, pCor):
    '''
    Arguments:
    * pLov, pCor * two profile objects
    Returns
    * Pres, Value, Qc * numpy ndarrays
    '''

    PresL, ValueL, QcL = pLov.read('CHLA',read_adjusted=True)
    if ('CHLA' in pCor.available_params) and  (pCor._my_float.status_var('CHLA')=='A'):
        PresC, ValueC, QcC = pCor.read('CHLA',read_adjusted=True)
    else:
        return PresL, ValueL, QcL
    if len(PresC)<5:

        if len(PresL)>5:
            print "few values in Coriolis: using LOV"
            return PresL, ValueL, QcL
        else:
            print "few values in either dataset for " + pLov._my_float.filename
            return None, None, None
    if has_drift(PresC, ValueC):
        print pCor._my_float.filename + " has drift"
        return None, None, None

    if len(ValueL)> 0:
        chl_lov_zero=ValueL[0]
    else:
        chl_lov_zero = None
    CHL = quenching(pCor, PresC, ValueC,chl_lov_zero)

    ii=(PresC >= 400) & (PresC <= 600)
    if ii.sum() > 0:
        shift = CHL[ii].mean()
        CHL = CHL - shift
    return PresC, CHL, QcC

def treating_coriolis(pCor):
    if pCor._my_float.status_var('CHLA')=='A':
        Pres,Value, Qc=pCor.read('CHLA', read_adjusted=True)
        if len(Pres)<5:
            print "few values in Coriolis for " + pCor._my_float.filename
            return None, None, None

        CHL = quenching(pCor, Pres, Value, None)
        ii=(Pres >= 400) & (Pres <= 600)
        if ii.sum() > 0:
            shift = CHL[ii].mean()
            CHL = CHL - shift
        return Pres, CHL, Qc
    else:
        print "R -- not dumped ", pCor._my_float.filename
        return None, None, None

if __name__=="__main__":

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

            if len(ValueC) < 5:
                print "few values in Coriolis for " + pCor._my_float.filename
                continue

            if len(ValueL)<  1:
                print "no LOV data"
                import sys
                sys.exit()

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

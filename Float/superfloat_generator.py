from instruments import bio_float
from instruments import lovbio_float
from commons.time_interval import TimeInterval
from basins.region import Rectangle
from commons import calculated_depths
import pylab as pl
import numpy as np
TI = TimeInterval('2015','2016','%Y')
R = Rectangle(-6,36,30,46)

def quenching(profile_obj, PresChl, Chl, chl_lov_zero):
    Quenched_Chl = Chl.copy()
    if profile_obj.time.month in range(1,13):
        PresT, Temp, _ = profile_obj.read('TEMP', read_adjusted=False)
        mld = calculated_depths.mld(Temp, PresT, zref=0)
        ii=PresChl<=mld
        
        nLev= ii.sum()+1
        quench_value1 = Chl[:nLev].max()
        nLev = ii.sum()
        quench_value2 = Chl[:nLev].max()
        diff1 = np.abs(quench_value1-chl_lov_zero)
        diff2 = np.abs(quench_value2-chl_lov_zero)
        if diff1<diff2:
            quench_value = quench_value1
            nLev =ii.sum()+1
        else:
            quench_value = quench_value2
            nLev = ii.sum()
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

#PROFILES_COR =   bio_float.FloatSelector('CHLA', TI, R)
PROFILES_LOV =lovbio_float.FloatSelector('CHLA', TI, R)

for ip, pLov in enumerate(PROFILES_LOV[80:]):
    pCor=bio_float.from_lov_profile(pLov)
    if pCor is not None:
        PresL, ValueL, QcL = pLov.read('CHLA',read_adjusted=True)
        PresC, ValueC, QcC = pCor.read('CHLA',read_adjusted=False)
        CHL = quenching(pCor, PresC, ValueC,ValueL[0])
        
        if not are_identical(CHL, ValueL):
        #    if not are_shifted(CHL, ValueL):
                print ip
                fig,ax =pl.subplots()
                ax.plot(ValueL,PresL,'r', label="LOV")
                ax.plot(CHL  , PresC,'b.', label="COR")
                ax.invert_yaxis()
                ax.grid()
                ax.legend()
                fig.show()
                #pl.savefig()
                import sys
                sys.exit()
         
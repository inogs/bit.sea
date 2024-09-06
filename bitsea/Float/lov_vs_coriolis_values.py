from instruments import lovbio_float
from instruments import bio_float
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import numpy as np
import matplotlib.pyplot as pl

TI = TimeInterval('2015','2016','%Y')
R = Rectangle(-6,36,30,46)

PROFILES_LOV=lovbio_float.FloatSelector('CHLA', TI, R)

for ip, pLov in enumerate(PROFILES_LOV[:1]):
    pCor = bio_float.from_lov_profile(pLov)
    if pCor is not None:
        PresL, ValueL, QcL = pLov.read('CHLA',read_adjusted=True)
        PresC, ValueC, QcC = pCor.read('CHLA',read_adjusted=False)
        fig,ax =pl.subplots()
        ax.plot(ValueL,PresL,'r', label="LOV")
        ax.plot(ValueC,PresC,'b.', label="COR")
        ax.invert_yaxis()
        ax.grid()
        ax.legend()
        fig.show()
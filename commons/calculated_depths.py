import numpy as np
def mld(Temperature, Pres, zref=10.0, deltaTemp=0.1):
    '''
    Definition of mld of Sprintall and Roemmich [1999], based on direct observation of more than 1000 profiles
    depth where deltaT = 0.1C
    '''
    index_pres=np.abs((Pres-zref)).argmin()
    absdiff_array = np.abs(Temperature - Temperature[index_pres])
    for k, absdiff in enumerate(absdiff_array):
        if absdiff > deltaTemp:
            break
    return Pres[k]

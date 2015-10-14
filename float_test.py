from basins.region import Region, Rectangle
from Float.float_manager import FloatSelector, TimeInterval

if __name__ == '__main__':
    var = 'NITRATE'
    TI = TimeInterval('20150520','20150830','%Y%m%d')
    R = Rectangle(-6,36,30,46)
    
    FLOAT_LIST=FloatSelector(var, TI, R)
    print FLOAT_LIST

    for TheFloat in FLOAT_LIST[:1]:
        PN,N = TheFloat.read(var)
        PS,S = TheFloat.read('PSAL')
        PT,T = TheFloat.read('TEMP')


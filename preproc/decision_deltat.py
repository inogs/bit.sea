import numpy as np
import os
mydtype=[('time','S8'), ('deltaT',np.float32),('K',np.int),('J',np.int),('I',np.int)]
INPUTDIR="DeltaT/"
cattedfile=INPUTDIR + "DeltaT.txt"
command="cat  " + INPUTDIR + "DeltaT_* > " + cattedfile 
os.system(command)


A=np.loadtxt(cattedfile, dtype=mydtype)

if A['deltaT'].min() > 600:
    print "namelist.init.600"
else:
    if A['deltaT'].min() > 450:
        print "namelist.init.450"
    else:
        print "namelist.init.300"
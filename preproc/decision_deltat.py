import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Prints in standard output the name of the choosen namelist.init file:
    One of these:
    - namelist.init.300
    - namelist.init.450
    - namelist.init.600
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = '''Directory where dailt DeltaT_*txt files are '''
                                )
    return parser.parse_args()


args = argument()
import numpy as np
import os
from commons.utils import addsep
mydtype=[('time','S8'), ('deltaT',np.float32),('K',np.int),('J',np.int),('I',np.int)]
INPUTDIR=addsep(args.inputdir)
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
        if A['deltaT'].min() > 360:
            print "namelist.init.360"
        else:
            print "namelist.init.300"


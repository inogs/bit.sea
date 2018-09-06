import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    ''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = None,
                                required = True,
                                help = "")
    return parser.parse_args()

args = argument()

import numpy as np
from profiler import MassimiliReader, Matchup_Manager, ALL_PROFILES, TL, BASEDIR, T_INT
from commons.mask import Mask
import basins.V2 as OGS
from instruments.var_conversions import MASSIMILIVARS as NUTRVARS
from commons.utils import addsep
from commons.layer import Layer

from commons.utils import writetable

M = Matchup_Manager(ALL_PROFILES, TL, BASEDIR)
N = MassimiliReader()

TheMask = Mask(args.maskfile)
nav_lev = TheMask.zlevels
OUTDIR = addsep(args.outdir)

LayerList = [Layer(0,1000), Layer(0,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300)]
SUBLIST=[OGS.nwm, OGS.tyr, OGS.lev, OGS.ion]

nSub, nLayer = len(SUBLIST), len(LayerList)
rows_names=[sub.name for sub in SUBLIST]
column_names = [layer.string() for layer in LayerList]


for modelvarname in ["N1p","N3n","P_l"]:
    fname = OUTDIR + modelvarname
    STAT = np.zeros((nSub,nLayer,3),np.float32)*np.nan
    NUMB = np.zeros((nSub,nLayer),np.int32)
    
    for isub, sub in enumerate(SUBLIST):
        Profilelist=N.Selector(NUTRVARS[modelvarname],T_INT,sub)
        nProfiles=len(Profilelist)    
        print sub.name, nProfiles
        Matchup_basin = M.getMatchups(Profilelist, nav_lev, modelvarname,read_adjusted=True)
        for ilayer, layer in enumerate(LayerList):
            m_layer = Matchup_basin.subset(layer)
            npoints = m_layer.number()
            NUMB[isub,ilayer] = npoints
            if npoints>3:
                print "npoints=", npoints
                STAT[isub,ilayer,0]=m_layer.bias()
                STAT[isub,ilayer,1]=m_layer.RMSE()
                STAT[isub,ilayer,2]=m_layer.correlation()

    writetable(fname + ".bias.txt", STAT[:,:,0], rows_names, column_names, fmt='%5.3f\t')
    writetable(fname + ".rmse.txt", STAT[:,:,1], rows_names, column_names, fmt='%5.3f\t')
    writetable(fname + ".corr.txt", STAT[:,:,2], rows_names, column_names, fmt='%5.3f\t')
    writetable(fname + ".numb.txt", NUMB       , rows_names, column_names, fmt='%1.0f\t')





        
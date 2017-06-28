#NIT-LAYER-Y-CLASS4-CLIM-BIAS
import numpy as np
from commons.layer import Layer
from matchup.statistics import matchup

INPUTDIR="../../static_clim/"
VARLIST=['N1p','N3n','O2o','Ac','DIC']
METRICvar = {'N1p':'PHO',
             'N3n':'NIT',
             'O2o':'DO',
             'Ac':'ALK',
             'DIC':'DIC'}



LayerList = [Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]

nLayers = len(LayerList)


for var in VARLIST:
    metric = METRICvar[var] + "-LAYER-Y-CLASS4-CLIM-"
    print ""
    print metric + "BIAS", metric + "RMSD"
    REF = np.load(INPUTDIR + var + "ref_clim.npy")
    MOD = np.load(INPUTDIR + var + "mod_clim.npy")
    for ilayer, layer in enumerate(LayerList):
        refsubs = REF[:,ilayer]
        modsubs = MOD[:,ilayer]
        bad = np.isnan(refsubs) | np.isnan(modsubs)
        m = matchup(modsubs[~bad], refsubs[~bad])
        print  m.bias(), m.RMSE()


from basins import V2
SUBlist = V2.Pred.basin_list

for var in ['Ac','DIC']:
    metric = METRICvar[var] + "-PROF-Y-CLASS4-CLIM-"
    print ""
    print metric + "BIAS", metric + "RMSD"
    REF = np.load(INPUTDIR + var + "ref_clim.npy")
    MOD = np.load(INPUTDIR + var + "mod_clim.npy")
    for isub, sub in enumerate(SUBlist):
        refsubs = REF[isub,:]
        modsubs = MOD[isub,:]
        bad = np.isnan(refsubs) | np.isnan(modsubs)
        m = matchup(modsubs[~bad], refsubs[~bad])
        print sub.name, m.bias(), m.RMSE()

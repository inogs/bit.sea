from commons import timerequestors
from basins import V2 as OGS
from instruments import bio_float
from instruments.var_conversions import FLOATVARS
from instruments.matchup_manager import Matchup_Manager
from commons.mask import Mask
from commons.layer import Layer
import numpy as np
from matchup.statistics import matchup


R=timerequestors.Weekly_req(2016,6,28)
TheMask = Mask('/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc')
BASEDIR="/gpfs/work/IscrC_MYMEDBIO/COPERNICUS/bit.sea/validation/online/NRT3/"
INPUTDIR="/pico/home/usera07ogs/a07ogs00/OPA/V2C/wrkdir/2/POSTPROC/AVE_FREQ_1/online_validation/PREVIOUS/TMP/"

LAYERLIST=[Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]
VARLIST = ['P_l','N3n','O2o']
nSub   = len(OGS.NRT3.basin_list)
nDepth = len(LAYERLIST)
nVar   = len(VARLIST)

BIAS    = np.zeros((nVar,nSub,nDepth), np.float32)
RMSE    = np.zeros((nVar,nSub,nDepth), np.float32)
NPOINTS = np.zeros((nVar,nSub,nDepth), np.float32)

TI =R.time_interval
M = Matchup_Manager(TI,INPUTDIR,BASEDIR)




for ivar, var in enumerate(VARLIST):
    for isub, sub in enumerate(OGS.NRT3):
        Profilelist = bio_float.FloatSelector(FLOATVARS[var], TI, sub)
        nProfiles = len(Profilelist)
        Matchup_object_list=[]
        for ip in range(nProfiles):
            floatmatchup =  M.getMatchups([Profilelist[ip]], TheMask.zlevels, var, read_adjusted=True)
            Matchup_object_list.append(floatmatchup)

        for ilayer, layer in enumerate(LAYERLIST):
            Matchup_object_list_layer=[]
            MODEL = []
            REF   = []
            for floatmatchup in Matchup_object_list:
                m_layer = floatmatchup.subset(layer)
                if m_layer.number() > 0:
                    Matchup_object_list_layer.append(m_layer)
                    REF.append(m_layer.Ref.mean())
                    MODEL.append(m_layer.Model.mean())

            NPOINTS[ivar,isub,ilayer] = len(MODEL)
            if len(MODEL) > 0:
                M_LAYER = matchup(np.array(MODEL), np.array(REF))
                BIAS[ivar,isub,ilayer] = M_LAYER.bias()
                RMSE[ivar,isub,ilayer] = M_LAYER.RMSE()
        
import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Calculates Chl statistics on matchups of Model and Sat.
    Main result are
    - BGC_CLASS4_CHL_RMS_SURF_BASIN
    - BGC_CLASS4_CHL_BIAS_SURF_BASIN
    of deliverable CMEMS-Med-biogeochemistry-ScCP-1.0.pdf

    Other similar results are also calculated and dumped in
    the a pickle file provided as argument.
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(   '--satdir', '-s',
                                type = str,
                                required =False,
                                default = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI/MONTHLY_V4/",
                                help = ''' Input satellite dir'''
                                )

    parser.add_argument(   '--moddir', '-m',
                                type = str,
                                required =False,
                                default = "/pico/scratch/userexternal/gbolzon0/RA_CARBO/RA_02/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/",
                                help = ''' Input model dir'''
                                )

    parser.add_argument(   '--outfile', '-o',
                                type = str,
                                required = True,
                                default = 'export_data_ScMYValidation_plan.pkl',
                                help = 'Output pickle file')

    parser.add_argument(   '--levsel', '-l',
                                type = float,
                                required = True,
                                default = 0.,
                                help = 'Level of mask for data selection')

    parser.add_argument(   '--dep', '-d',
                                type = float,
                                required = False,
                                default = 10.,
                                help = 'Level of bottom layer for average')


    return parser.parse_args()


args = argument()



from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import commons.IOnames as IOnames
import numpy as np
import Sat.SatManager as Sat
import matchup.matchup as matchup
from commons.dataextractor import DataExtractor
from layer_integral.mapbuilder import MapBuilder
from commons.mask import Mask
from commons.submask import SubMask
from basins import OGS
from commons.layer import Layer
import pickle
import os

def weighted_mean(Conc, Weight):

    Weight_sum      = Weight.sum()
    Mass            = (Conc * Weight).sum()
    Weighted_Mean   = Mass/Weight_sum
    return Weighted_Mean

maskfile = os.getenv('MASKFILE')
TheMask = Mask(maskfile)

Timestart=os.getenv("START_DATE")
Time__end=os.getenv("END_DATE")

MODEL_DIR= args.moddir
REF_DIR  = args.satdir
outfile  = args.outfile


TI    = TimeInterval(Timestart,Time__end,"%Y%m%d")

sat_TL   = TimeList.fromfilenames(TI, REF_DIR  ,"*.nc", prefix="", dateformat="%Y%m%d")
model_TL = TimeList.fromfilenames(TI, MODEL_DIR,"*.nc")

IOname = IOnames.IOnames('IOnames_sat.xml')

nFrames = model_TL.nTimes
nSUB = len(OGS.P.basin_list)

jpk,jpj,jpi =TheMask.shape
dtype = [(sub.name, np.bool) for sub in OGS.P]
SUB = np.zeros((jpj,jpi),dtype=dtype)
for sub in OGS.P:
    sbmask         = SubMask(sub,maskobject=TheMask).mask
    SUB[sub.name]  = sbmask[0,:,:]

masksel_2D = TheMask.mask_at_level(args.levsel)

BGC_CLASS4_CHL_RMS_SURF_BASIN      = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN     = np.zeros((nFrames,nSUB),np.float32)
MODEL_MEAN                         = np.zeros((nFrames,nSUB),np.float32)
SAT___MEAN                         = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG  = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG = np.zeros((nFrames,nSUB),np.float32)

# This is the surface layer choosen to match satellite chl data
surf_layer = Layer(0,args.dep)

for itime, modeltime in enumerate(model_TL.Timelist):
    print modeltime
    CoupledList = sat_TL.couple_with([modeltime])
    sattime = CoupledList[0][0]
    satfile = REF_DIR + sattime.strftime(IOname.Input.dateformat) + IOname.Output.suffix + ".nc"
    modfile = model_TL.filelist[itime]

    De         = DataExtractor(TheMask,filename=modfile, varname='P_i')
    Model      = MapBuilder.get_layer_average(De, surf_layer)
    #ncIN = NC.netcdf_file(modfile,'r')
    #Model = ncIN.variables['P_i'].data[0,0,:,:].copy()#.astype(np.float64)
    #Model = ncIN.variables['lchlm'].data.copy()
    #ncIN.close()

    Sat16 = Sat.readfromfile(satfile,var='lchlm') #.astype(np.float64)


    cloudsLand = (np.isnan(Sat16)) | (Sat16 > 1.e19) | (Sat16<0)
    modelLand  = np.isnan(Model) #lands are nan
    nodata     = cloudsLand | modelLand
    selection = ~nodata & masksel_2D
    M = matchup.matchup(Model[selection], Sat16[selection])

    for isub, sub in enumerate(OGS.P):
        selection = SUB[sub.name] & (~nodata) & masksel_2D
        M = matchup.matchup(Model[selection], Sat16[selection])
        BGC_CLASS4_CHL_RMS_SURF_BASIN[itime,isub]  = M.RMSE()
        BGC_CLASS4_CHL_BIAS_SURF_BASIN[itime,isub] = M.bias()
        weight = TheMask.area[selection]
        MODEL_MEAN[itime,isub] = weighted_mean( M.Model,weight)
        SAT___MEAN[itime,isub] = weighted_mean( M.Ref,  weight)

        Mlog = matchup.matchup(np.log10(Model[selection]), np.log10(Sat16[selection])) #add matchup based on logarithm
        BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[itime,isub]  = Mlog.RMSE()
        BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[itime,isub] = Mlog.bias()

BGC_CLASS4_CHL_EAN_RMS_SURF_BASIN  = BGC_CLASS4_CHL_RMS_SURF_BASIN.mean(axis=0)
BGC_CLASS4_CHL_EAN_BIAS_SURF_BASIN = BGC_CLASS4_CHL_BIAS_SURF_BASIN.mean(axis=0)


LIST   =[i for i in range(7)]

LIST[0]=model_TL.Timelist
LIST[1]=BGC_CLASS4_CHL_RMS_SURF_BASIN
LIST[2]=BGC_CLASS4_CHL_BIAS_SURF_BASIN
LIST[3]=MODEL_MEAN
LIST[4]=SAT___MEAN
LIST[5]=BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG
LIST[6]=BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG

fid = open(outfile,'wb')
pickle.dump(LIST, fid)
fid.close()

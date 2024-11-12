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
    parser.add_argument(   '--inputmodeldir', '-i',
                                type = str,
                                required =True,
                                help = ''' Input model dir, where P_l files are, usually ../wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/'''
                                )
    parser.add_argument(   '--outfile', '-o',
                                type = str,
                                required = True,
                                default = 'export_data_ScMYValidation_plan.pkl',
                                help = 'Output pickle file')
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')
    parser.add_argument(   '--coastness', '-c',
                                type = str,
                                required = True,
                                choices = ['coast','open_sea','everywhere'],
                                help = 'definition of mask to apply to the statistics')



    return parser.parse_args()


args = argument()



from bitsea.commons.Timelist import TimeList
from bitsea.commons.time_interval import TimeInterval
import numpy as np
import os
import bitsea.Sat.SatManager as Sat
import bitsea.matchup.matchup as matchup
from bitsea.commons.dataextractor import DataExtractor
from bitsea.layer_integral.mapbuilder import MapBuilder
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
from bitsea.basins import V2 as OGS
from bitsea.commons.layer import Layer
from bitsea.commons.utils import addsep
import pickle

def weighted_mean(Conc, Weight):

    Weight_sum      = Weight.sum()
    Mass            = (Conc * Weight).sum()
    Weighted_Mean   = Mass/Weight_sum
    return Weighted_Mean

TheMask = Mask.from_file(args.maskfile)
Sup_mask = TheMask.cut_at_level(0)
MODEL_DIR= addsep(args.inputmodeldir)
REF_DIR  = addsep(args.satdir)
outfile  = args.outfile


Timestart="20141205"
Time__end="20500901"
TI    = TimeInterval(Timestart,Time__end,"%Y%m%d")
dateformat ="%Y%m%d"
print(REF_DIR)
sat_TL   = TimeList.fromfilenames(TI, REF_DIR  ,"*.nc", prefix="", dateformat=dateformat)
print(sat_TL.nTimes)
model_TL = TimeList.fromfilenames(TI, MODEL_DIR,"ave.*00:00*P_l.nc")
suffix = os.path.basename(sat_TL.filelist[0])[8:]


nFrames = model_TL.nTimes
nSUB = len(OGS.P.basin_list)

jpk,jpj,jpi =TheMask.shape
dtype = [(sub.name, bool) for sub in OGS.P]
SUB = np.zeros((jpj,jpi),dtype=dtype)
for sub in OGS.Pred:
    print sub.name
    sbmask         = SubMask(sub,maskobject=Sup_mask).mask
    SUB[sub.name]  = sbmask[0,:,:]

for sub in OGS.Pred:
    SUB['med'] = SUB['med'] | SUB[sub.name]

mask200_2D = TheMask.mask_at_level(200.0)
mask0_2D = TheMask.mask_at_level(0.0)
if args.coastness == 'coast':
    coastmask=mask0_2D & (~mask200_2D)
if args.coastness == "open_sea"  : coastmask = mask200_2D
if args.coastness == "everywhere": coastmask = mask0_2D

BGC_CLASS4_CHL_RMS_SURF_BASIN      = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN     = np.zeros((nFrames,nSUB),np.float32)
MODEL_MEAN                         = np.zeros((nFrames,nSUB),np.float32)
SAT___MEAN                         = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG  = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG = np.zeros((nFrames,nSUB),np.float32)
MODEL_MEAN_allpoints               = np.zeros((nFrames,nSUB),np.float32)

# This is the surface layer choosen to match satellite chl data
surf_layer = Layer(0,10)

for itime, modeltime in enumerate(model_TL.Timelist):
    print modeltime
    CoupledList = sat_TL.couple_with([modeltime])
    print(CoupledList)
    sattime = CoupledList[0][0]
    satfile = REF_DIR + sattime.strftime(dateformat) + suffix
    modfile = model_TL.filelist[itime]

    De         = DataExtractor(TheMask,filename=modfile, varname='P_l')
    Model      = MapBuilder.get_layer_average(De, surf_layer)
    #ncIN = NC.netcdf_file(modfile,'r')
    #Model = ncIN.variables['P_i'].data[0,0,:,:].copy()#.astype(np.float64)
    #Model = ncIN.variables['lchlm'].data.copy()
    #ncIN.close()

    Sat16 = Sat.readfromfile(satfile,var='CHL') #.astype(np.float64)


    cloudsLand = (np.isnan(Sat16)) | (Sat16 > 1.e19) | (Sat16<0)
    modelLand  = np.isnan(Model) #lands are nan
    nodata     = cloudsLand | modelLand
    selection = ~nodata & coastmask
    M = matchup.matchup(Model[selection], Sat16[selection])

    for isub, sub in enumerate(OGS.P):
        selection = SUB[sub.name] & (~nodata) & coastmask
        M = matchup.matchup(Model[selection], Sat16[selection])
        BGC_CLASS4_CHL_RMS_SURF_BASIN[itime,isub]  = M.RMSE()
        BGC_CLASS4_CHL_BIAS_SURF_BASIN[itime,isub] = M.bias()
        weight = TheMask.area[selection]
        MODEL_MEAN[itime,isub] = weighted_mean( M.Model,weight)
        SAT___MEAN[itime,isub] = weighted_mean( M.Ref,  weight)

        Mlog = matchup.matchup(np.log10(Model[selection]), np.log10(Sat16[selection])) #add matchup based on logarithm
        BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[itime,isub]  = Mlog.RMSE()
        BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[itime,isub] = Mlog.bias()
        selectionall = SUB[sub.name] & coastmask
        MODEL_MEAN_allpoints[itime,isub] = np.nanmean(Model[selectionall])

BGC_CLASS4_CHL_EAN_RMS_SURF_BASIN  = BGC_CLASS4_CHL_RMS_SURF_BASIN.mean(axis=0)
BGC_CLASS4_CHL_EAN_BIAS_SURF_BASIN = BGC_CLASS4_CHL_BIAS_SURF_BASIN.mean(axis=0)


LIST   =[i for i in range(8)]

LIST[0]=model_TL.Timelist
LIST[1]=BGC_CLASS4_CHL_RMS_SURF_BASIN
LIST[2]=BGC_CLASS4_CHL_BIAS_SURF_BASIN
LIST[3]=MODEL_MEAN
LIST[4]=SAT___MEAN
LIST[5]=BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG
LIST[6]=BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG
LIST[7]=MODEL_MEAN_allpoints

fid = open(outfile,'wb')
pickle.dump(LIST, fid)
fid.close()

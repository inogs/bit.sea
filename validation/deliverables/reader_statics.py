import pickle
import numpy as np

fid = open('export_data_ScMYValidation_plan_statics.pkl')
LIST = pickle.load(fid)
fid.close()

BGC_CLASS4_PHOS_RMS_LAYER_BASIN = LIST[0]
BGC_CLASS4_NIT_RMS_LAYER_BASIN  = LIST[1]
BGC_CLASS4_O2_RMS_LAYER_BASIN   = LIST[2]
SUBlist                         = LIST[3]
LAYERLIST                       = LIST[4]
#QUID TABLE IV.4
np.savetxt('RMS_PHOS.csv',LIST[0],fmt='%10.5f',delimiter=',')
np.savetxt('RMS__NIT.csv',LIST[1],fmt='%10.5f',delimiter=',')
np.savetxt('RMS___O2.csv',LIST[2],fmt='%10.5f',delimiter=',')



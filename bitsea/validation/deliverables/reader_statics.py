import pickle
import numpy as np
import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Plot 4.3 and 4.4
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required = True,
                            default = './table4.3/',
                            help = 'Out directory')


    return parser.parse_args()

args = argument()

fid = open('export_data_ScMYValidation_plan_statics.pkl')
LIST = pickle.load(fid)
fid.close()

BGC_CLASS4_PHOS_RMS_LAYER_BASIN = LIST[0]
BGC_CLASS4_NIT_RMS_LAYER_BASIN  = LIST[1]
BGC_CLASS4_O2_RMS_LAYER_BASIN   = LIST[2]
SUBlist                         = LIST[3]
LAYERLIST                       = LIST[4]
#QUID TABLE IV.4
np.savetxt(args.outdir+'/RMS_PHOS.csv',LIST[0],fmt='%10.5f',delimiter=',')
np.savetxt(args.outdir+'/RMS__NIT.csv',LIST[1],fmt='%10.5f',delimiter=',')
np.savetxt(args.outdir+'/RMS___O2.csv',LIST[2],fmt='%10.5f',delimiter=',')

import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    requires frequency of float DA
    ''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(   '--template',"-t",
                                type = str,
                                required = True,
                                help = 'path of the template xml file')
    parser.add_argument(   '--ref_daTimes','-r',
                                type = str,
                                required = True,
                                help = 'output or merge_DaTimes')
    parser.add_argument(   '--daTimes_sat','-s',
                                type = str,
                                required = True,
                                help = 'output or merge_DaTimes')
    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                help = '/some/path or wrkdir/DA__FREQ_1')
    parser.add_argument(   '--daTimes_dir','-d',
                                type = str,
                                required = True,
                                help = 'output dir of profiles_dates_DAfreq.py')
    return parser.parse_args()

args = argument()

from commons.utils import file2stringlist
import numpy as np
from commons.utils import addsep

def dump_template(ORIG, outfile,SAT_OBS,ARGO,ASS_P_l,ASS_N3n,ASS_O2o):
    LINES=[]
    for line in ORIG:
        newline=line
        if (line.find("@@SAT_OBS@@") != -1):  newline=line.replace("@@SAT_OBS@@",SAT_OBS)
        if (line.find("@@ARGO@@")   != -1):   newline=line.replace("@@ARGO@@"   ,  ARGO)
        if (line.find("@@ASS_P_l@@") != -1):  newline=line.replace("@@ASS_P_l@@",np.str(np.int(ASS_P_l)+np.int(SAT_OBS)))
        if (line.find("@@ASS_N3n@@") != -1):  newline=line.replace("@@ASS_N3n@@",ASS_N3n)
        if (line.find("@@ASS_O2o@@") != -1):  newline=line.replace("@@ASS_O2o@@",ASS_O2o)
        ASS_P_lN3n = '0'
        if (ASS_P_l=='1') and (ASS_N3n=='0') :
            ASS_P_lN3n = '1'
        if (line.find("@@ASS_P_lN3n@@") != -1):  newline=line.replace("@@ASS_P_lN3n@@",np.str(np.int(ASS_P_lN3n)))
        ASS_N3nO2o = '0'
        if (ASS_N3n=='1') or (ASS_O2o=='1') :
            ASS_N3nO2o = '1'
        if (line.find("@@ASS_N3nO2o@@") != -1):  newline=line.replace("@@ASS_N3nO2o@@",np.str(np.int(ASS_N3nO2o)))
        LINES.append(newline + "\n")
    fid=open(filename,"w")
    fid.writelines(LINES)
    fid.close()

DIR   =addsep(args.daTimes_dir)
OUTDIR=addsep(args.outdir)

ORIG          = file2stringlist(args.template)
dateDA        = file2stringlist(args.ref_daTimes)
dateDAsat     = file2stringlist(args.daTimes_sat)

dateDAfloat    = file2stringlist(DIR + "daTimes_float_N3nP_lO2o")
dateDAfloatP_l = file2stringlist(DIR + "daTimes_float_P_l")
dateDAfloatN3n = file2stringlist(DIR + "daTimes_float_N3n")
dateDAfloatO2o = file2stringlist(DIR + "daTimes_float_O2o")


for d in dateDA:
    filename=OUTDIR + "satfloat." + d[:] + ".nml"
    ARGO    = str(int(d in dateDAfloat)) # "1" if True, "0" if False
    SAT_OBS = str(int(d in dateDAsat))
    ASS_P_l = str(int(d in dateDAfloatP_l))
    ASS_N3n = str(int(d in dateDAfloatN3n))
    ASS_O2o = str(int(d in dateDAfloatO2o))
    dump_template(ORIG, filename, SAT_OBS, ARGO, ASS_P_l, ASS_N3n, ASS_O2o)


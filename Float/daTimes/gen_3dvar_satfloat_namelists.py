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

def dump_template(ORIG, outfile,SAT_OBS,ARGO,ASS_P_l,ASS_N3n):
    LINES=[]
    for line in ORIG:
        newline=line
        if (line.find("@@SAT_OBS@@") != -1):  newline=line.replace("@@SAT_OBS@@",SAT_OBS)
        if (line.find("@@ARGO@@")   != -1):   newline=line.replace("@@ARGO@@"   ,  ARGO)
        if (line.find("@@ASS_P_l@@") != -1):  newline=line.replace("@@ASS_P_l@@",np.str(np.int(ASS_P_l)+np.int(SAT_OBS)))
        if (line.find("@@ASS_N3n@@") != -1):  newline=line.replace("@@ASS_N3n@@",ASS_N3n)
        if (line.find("@@ASS_P_lN3n@@") != -1):  newline=line.replace("@@ASS_P_lN3n@@",np.str(np.int(ASS_P_l)-np.int(ASS_N3n)))
        LINES.append(newline + "\n")
    fid=open(filename,"w")
    fid.writelines(LINES)
    fid.close()

DIR   =addsep(args.daTimes_dir)
OUTDIR=addsep(args.outdir)

ORIG          = file2stringlist(args.template)
dateDA        = file2stringlist(args.ref_daTimes)
dateDAsat     = file2stringlist(args.daTimes_sat)

dateDAfloat    = file2stringlist(DIR + "daTimes_float_N3nP_l")
dateDAfloatP_l = file2stringlist(DIR + "daTimes_float_P_l")
dateDAfloatN3n = file2stringlist(DIR + "daTimes_float_N3n")


for d in dateDA:
    filename=OUTDIR + "satfloat." + d[:] + ".nml"
    ARGO    = str(int(d in dateDAfloat)) # "1" if True, "0" if False
    SAT_OBS = str(int(d in dateDAsat))
    ASS_P_l = str(int(d in dateDAfloatP_l))
    ASS_N3n = str(int(d in dateDAfloatN3n))
    dump_template(ORIG, filename, SAT_OBS, ARGO, ASS_P_l, ASS_N3n)



from commons.utils import file2stringlist
import numpy as np

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



ORIG=file2stringlist('satfloat.template')
dateDAfloat = file2stringlist("daFloat")
dateDAsat   = file2stringlist("daSat")
dateDA      = file2stringlist("daTimes")
dateDAfloatP_l = file2stringlist("daFloat_P_l")
dateDAfloatN3n = file2stringlist("daFloat_N3n")

OUTDIR="OUT/"
for d in dateDA:
    filename=OUTDIR + "satfloat." + d[:] + ".nml"
    #print filename
    ARGO    = str(int(d in dateDAfloat)) # "1" if True, "0" if False
    SAT_OBS = str(int(d in dateDAsat))
    ASS_P_l = str(int(d in dateDAfloatP_l))
    ASS_N3n = str(int(d in dateDAfloatN3n))
    dump_template(ORIG, filename, SAT_OBS, ARGO, ASS_P_l, ASS_N3n)






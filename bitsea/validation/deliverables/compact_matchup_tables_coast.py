import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    plot something
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            help = ''' Output dir'''
                            )

    parser.add_argument(   '--indir', '-i',
                            type = str,
                            required =True,
                            help = ''' Input dir'''
                            )

    parser.add_argument(   '--variable', '-v',
                            type = str,
                            required = True,
                            help = 'variable name')

    return parser.parse_args()

args = argument()

import numpy as np
from commons.utils import addsep
from commons.utils import writetable

INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)
var = args.variable

RMSE = np.loadtxt(INDIR + '/' + var + '.rmse.txt',
            skiprows=1,usecols=[1,2])

subnames = np.loadtxt(INDIR + '/' + var + '.rmse.txt',
            dtype='S',usecols=0)

header = [
    '0-60m', 
    '60-200m',
    'nprof/nobs',
    ]

metrics = np.zeros((RMSE.shape[0],len(header)-1))
# metrics = np.zeros((RMSE.shape[0],len(header)))
metrics[:,:] = RMSE

NPRF = np.loadtxt(INDIR + '/' + var + '.nprofiles.txt',
            skiprows=1,usecols=1)

NOBS = np.loadtxt(INDIR + '/' + var + '.numb.txt',
            skiprows=1,usecols=3)

np_nobs = np.zeros((len(subnames),1),dtype='S10')
for jj in range(len(subnames)):
    np_nobs[jj] = '%d' %NPRF[jj] + '/' + '%d' %NOBS[jj]

M = np.concatenate((metrics,np_nobs),axis=1)
filename = OUTDIR + '/coast.' + var + '.txt'
writetable(filename, M, subnames, header, fmt='%5.3f\t'*2 + '%s\t')
# writetable(filename, metrics, subnames, header)

# filename = OUTDIR + '/coast_np_nobs.' + var + '.txt'
# np.savetxt(filename,np_nobs,fmt='%s')
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
            skiprows=1,usecols=[ii for ii in range(1,8)])

CORR = np.loadtxt(INDIR + '/' + var + '.corr.txt',
            skiprows=1,usecols=8)

subnames = np.loadtxt(INDIR + '/' + var + '.corr.txt',
            dtype='S',usecols=0)

header = [
    '0-30m', 
    '30-60m',
    '60-100m',
    '100-150m',
    '150-300m',
    '300-600m',
    '600-1000m',
    'CORR',
    'nprof/nobs'
    ]

metrics = np.zeros((RMSE.shape[0],len(header)-1))
# metrics = np.zeros((RMSE.shape[0],len(header)))
metrics[:,0:7] = RMSE
metrics[:,7] = CORR

NPRF = np.loadtxt(INDIR + '/' + var + '.nprofiles.txt',
            skiprows=1,usecols=1)

NOBS = np.loadtxt(INDIR + '/' + var + '.numb.txt',
            skiprows=1,usecols=8)

np_nobs = np.zeros((len(subnames),1),dtype='S10')
for jj in range(len(subnames)):
    np_nobs[jj] = '%d' %NPRF[jj] + '/' + '%d' %NOBS[jj]

M = np.concatenate((metrics,np_nobs),axis=1)
filename = OUTDIR + '/open.' + var + '.txt'
writetable(filename, M, subnames, header, fmt='%5.3f\t'*8 + '%s\t')
# writetable(filename, metrics, subnames, header)

# filename = OUTDIR + '/open_np_nobs.' + var + '.txt'
# np.savetxt(filename,np_nobs,fmt='%s')
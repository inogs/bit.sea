import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Calculate eof for Mediterranean open sea and coastal sub-basins
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Outputdir'''
                            )

    parser.add_argument(   '--instatp', '-i',
                            type = str,
                            required = True,
                            help = 'Input STAT_PROFILES'
                            )

    parser.add_argument(   '--maskfile', '-m',
                            type = str,
                            required = True,
                            help = 'model maskfile meshmask.nc'
                            )

    parser.add_argument(   '--varname', '-v',
                            type = str,
                            required = True,
                            help = 'varname '
                            )

    parser.add_argument(   '--eofdepth', '-d',
                            type = float,
                            required = True,
                            help = 'EOF depth [m] '
                            )


    return parser.parse_args()

args = argument()

import numpy as np
import pickle as pkl
from eofs.standard import Eof
from commons.utils import addsep
from commons.timerequestors import Clim_month
from commons.mask import Mask
from basins import V0_adr as OGS



INDIR = addsep(args.instatp)
OUTDIR = addsep(args.outdir)
varname = args.varname
maskfile = args.maskfile
eofdepth = args.eofdepth

TheMask = Mask(maskfile)
zlevs = TheMask.zlevels
kdep = TheMask.getDepthIndex(eofdepth)+1

nSub = len(OGS.Pred.basin_list)
neof = kdep+1


doread = True
if doread:
    print 'do read'
    filepkl = INDIR + varname + '.pkl'
    fid = open(filepkl,'r')
    (profiles,dates) = pkl.load(fid)
    fid.close()


if eofdepth<=200:
    DICTind = {
        'open': 1,
        'coast1': 3,
        'coast2': 4,
    }
else:
    DICTind = {
        'open': 1,
    }

for im in range(1,13):
    print im
    req = Clim_month(im)
    ii,_ = dates.select(req)
    nFrames = len(ii)
    if nFrames<neof:
        print 'n dates < of n levels'
        print 'EXIT beacuse the script not supposed to work'
        import sys
        sys.exit(0)
    for aType in DICTind.keys():
        print ' ... ' + aType
        evc = np.zeros((nSub,neof,kdep+1))
        eva = np.zeros((nSub,neof))
        mprof = profiles[ii,:,DICTind[aType],:kdep+1,0]
        for kk in range(1,kdep+1):
            if np.any(np.isnan(mprof[:,:,kk])):
                nans = np.isnan(mprof[:,:,kk])
                mprof[nans,kk] = mprof[nans,kk-1]
        if np.any(np.isnan(mprof)):
            print 'Some nan in profiles'
            import sys
            sys.exit(0)
        meanP = np.nanmean(mprof,0)
        meanPmat = np.zeros_like(mprof)
        for iframe in range(nFrames):
            meanPmat[iframe,:,:] = meanP
        anomalies = mprof-meanPmat 
        for isub in range(nSub):
            anosub = anomalies[:,isub,:]
            sol = Eof(anosub)
            evc[isub,:,:] = sol.eofs()[:neof,:kdep+1]
            eva[isub,:] = (sol.eigenvalues()[:neof])**.5

        filepkl = OUTDIR + '/evc' + aType + '.%02d.pkl' %(im)  
        fid = open(filepkl,'wb')
        pkl.dump(evc,fid)
        fid.close()

        filepkl = OUTDIR + '/eva' + aType + '.%02d.pkl' %(im)  
        fid = open(filepkl,'wb')
        pkl.dump(eva,fid)
        fid.close()



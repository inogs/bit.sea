import argparse

def argument():
    parser = argparse.ArgumentParser(description='''
    complete weekly maps using monthly mean or climatology at model resolution
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--maskfile', '-m',
                        type=str,
                        default=None,
                        required=True,
                        help=''' Path of maskfile''')

    parser.add_argument('--indir', '-i',
                        type=str,
                        default=None,
                        required=True,
                        help=" Path of gap-filled weekly files")

    parser.add_argument('--dirind', '-n',
                        type=str,
                        default=None,
                        required=True,
                        help=" Path of indexes files")

    parser.add_argument('--perc', '-p',
                        type=str,
                        default=None,
                        required=True,
                        help=" Perc to increas kd")

    parser.add_argument('--outdir', '-o',
                        type=str,
                        default=None,
                        required=True,
                        help=" outdir")

    return parser.parse_args()

args = argument()


import numpy as np
import scipy.io.netcdf as NC
import glob,os
from eofs.standard import Eof
from commons.mask import Mask
from commons.utils import addsep

def modulation_func(myx,mycorr):
    mycenter=0.2*0.5
    return (mycorr-1.0) * (0.5 * np.pi - np.arctan(100.0*(myx-mycenter)))/np.pi + 1.0


INPDIR = addsep(args.indir)
DIRIND = addsep(args.dirind)
OUTDIR = addsep(args.outdir)
perc = np.float(args.perc)
TheMask = Mask(args.maskfile)
Lon = TheMask.xlevels
Lat = TheMask.ylevels

tmask = TheMask.mask

_,jpj,jpi = TheMask.shape

tmaskS = tmask[0,:]
nump = np.sum(tmaskS)


inds,_ = TheMask.convert_lon_lat_to_indices(-5,36)
inde,_ = TheMask.convert_lon_lat_to_indices(-2,36)


filelist = glob.glob(INPDIR + '*.npy')
filelist.sort()
nsample = len(filelist)

matrixMaps = np.zeros((nsample,nump))

print('Loading data and evaluating anonalies')
ii = 0
for infile in filelist:
    datec = os.path.basename(infile)[0:8]
    map = np.load(infile)
    vec = map[map>0]
    matrixMaps[ii,:] = vec
    ii = ii+1

meanMaps = np.mean(matrixMaps,0)

matrixAnom = np.zeros((nsample,nump))
for ii in range(nsample):
    matrixAnom[ii,:] = matrixMaps[ii,:]-meanMaps

#matrixAnom = recLim0
print('Eof evaluation and field reconstruction')
ans = Eof(matrixAnom)

eofm = ans.eofs()
eigm = ans.eigenvalues()

exvar = 0
found = False
limvar = 99
for ii in range(nsample):
    exvar = ans.varianceFraction()[ii]*100+exvar
    if (exvar>limvar) & ~found:
       neofLim = ii+1
       found = True
       varii = exvar



recAno = ans.reconstructedField(neofs=neofLim)

def write_nc(outfile,var2save,var):
    ncOUT   = NC.netcdf_file(outfile,"w")

    ncOUT.createDimension('x',jpi)
    ncOUT.createDimension('y',jpj)

    ncvar = ncOUT.createVariable(var,'f',('y','x'))
    ncvar[:] = var2save
    ncvar = ncOUT.createVariable('nav_lon','f',('y','x'))
    ncvar[:] = Lon
    ncvar = ncOUT.createVariable('nav_lat','f',('y','x'))
    ncvar[:] = Lat

    setattr(ncOUT.variables[var],'missing_value',1.e+20);

    ncOUT.close()

print('Writing NC files')

indI = np.load(DIRIND + '/indI.npy')
indJ = np.load(DIRIND + '/indJ.npy')
var = 'kextfact'
factinc = 1 + perc/100.
for id in range(nsample):
    filenpy = filelist[id]
    datec = os.path.basename(filenpy)[0:8]
    print(datec)
    outfile = OUTDIR + '/KextF_' + datec + '-00:00:00.nc'
    var2save = np.zeros((jpj,jpi))
    for ip in range(nump):
        ii = indI[ip]
        jj = indJ[ip]
        recFld = recAno[id,ip]+meanMaps[ip]
        var2save[jj,ii] = recFld *  modulation_func(recFld,factinc)
#       var2save[jj,ii] = recFld * factinc
    maskv = var2save==0
    var2save[maskv] = np.nan
    meanAlb = np.nanmean(var2save[:,inds:inde])
    var2save[:,:inds] = meanAlb
    var2save[maskv] = 1.e+20
    write_nc(outfile,var2save,var)
    



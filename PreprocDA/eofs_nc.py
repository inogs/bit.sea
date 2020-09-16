import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Eofs calibration with variance for open sea and coast
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Outputdir'''
                            )

    parser.add_argument(   '--ineof', '-i',
                            type = str,
                            required = True,
                            help = 'Coastal and open sea eofs pkl dir'
                            )

    parser.add_argument(   '--maskfile', '-m',
                            type = str,
                            required = True,
                            help = 'model maskfile meshmask.nc'
                            )

    parser.add_argument(   '--coastnpy', '-c',
                            type = str,
                            required = True,
                            help = 'npy coast k-regions mask'
                            )

    parser.add_argument(   '--grid3dvar', '-g',
                            type = str,
                            required = True,
                            help = '3dvar grid'
                            )

    parser.add_argument(   '--vardir', '-v',
                            type = str,
                            required = True,
                            help = 'varnod dir'
                            )

    parser.add_argument(   '--neof', '-n',
                            type = str,
                            required = True,
                            help = 'varnod dir'
                            )

    return parser.parse_args()

args = argument()



import numpy as np
from commons.utils import addsep
from commons.mask import Mask
from commons.submask import SubMask
from basins import V0_adr as OGS
import pickle as pkl
import scipy.io.netcdf as NC
#from maskload import tmask,SUB,SUBlist,jpk,jpj,jpi

CASES = ['ext0','ext1','prfc']

EOFDIR = addsep(args.ineof)
VARDIR = addsep(args.vardir)
OUTEOF = addsep(args.outdir)
mapcfile = args.coastnpy
maskfile = args.maskfile
var3dgridfile = args.grid3dvar
neof = int(args.neof)

TheMask = Mask(maskfile)
ind200 = TheMask.getDepthIndex(200)+1
jpk,jpj,jpi = TheMask.shape
mask200 = TheMask.mask_at_level(200)
masks = TheMask.mask_at_level(0)
mask_open = masks & mask200
maskcoast = masks & (mask200==False)


mapc = np.load(mapcfile)
if (mapc.shape==(jpj,jpi))==False:
    print 'incorrect maskcoast shape EXIT'
    import sys
    sys.exit(0)
#version = 'v02' # version of model error limit
                # v01 : 1/2
                # v02 : 3/4
#GROUPTYPE = 'MONTHGRCLIM'
#MAPCDIR = '../../OUTPUT/' + GROUPTYPE + '/CLASS_SERIES/'
#MAPCDIR = 'MAPCinterp24/'
#EOFCDIR = '../../OUTPUT/' + GROUPTYPE + '/EOFCOAST/'
#EVADIR = '../../OUTPUT/CMEMScoast/' + GROUPTYPE + '/EOFCOAST600/'



MM = NC.netcdf_file(var3dgridfile,'r')
nregs = int(np.max(MM.variables['regs'].data))
regs = MM.variables['regs'].data.copy()
tgr = MM.variables['tmsk'].data[0,:,:].copy()

def write_eof(outfile,eofp,eofa):
    ncOUT   = NC.netcdf_file(outfile,"w")

    ncOUT.createDimension('nreg',nregs)
    ncOUT.createDimension('nlev',jpk)
    ncOUT.createDimension('neof',neof)

    ncvar = ncOUT.createVariable('eva','f',('neof','nreg'))
    ncvar[:] = eofa;
    ncvar = ncOUT.createVariable('evc','f',('neof','nlev','nreg'))
    ncvar[:] = eofp;

    ncOUT.close()



AREASlist = [sub.name for sub in OGS.Pred]
dtype = [(sub.name, np.bool) for sub in OGS.Pred]
SUBareas = np.zeros((jpj,jpi),dtype=dtype)
SUBtot = np.zeros((jpj,jpi),np.int)
for isub,sub in enumerate(OGS.Pred):
    print sub.name
    SUBareas[sub.name] = SubMask(sub,maskobject=TheMask).mask_at_level(0)
    SUBtot[SUBareas[sub.name]] = isub
SUBtot[tgr==0] = -1000


mask2 = maskcoast &(SUBtot>=0)
[J2,I2] = np.where(mask2)
nr2 = I2.size

#mask0 = mask_open & (SUBtot)
mask0 = masks & (SUBtot>=0)# mapc interpolated not all point definedin the meshmask
[J0,I0] = np.where(mask0)
nr0 = I0.size


for im in range(12):
    isopen = np.zeros(nregs,dtype=bool)
    mstr = '%02d' %(im+1)
    print(mstr)

#    EOF = NC.netcdf_file(EOFDIR + 'eof.9.'+ mstr + '.nc','r')
    VAR = NC.netcdf_file(VARDIR + 'var2D.' + mstr + '.nc')
    varmap = VAR.variables['variance'].data[0,:,:].copy()
    VAR.close()
    #mapc[mapc>1.e+19] = np.nan

    eofp = np.zeros((neof,jpk,nregs))
    eofa = np.zeros((neof,nregs))

    evc = {}
    eva = {}
    for aType in ['open','coast1','coast2']:
        fileevc = EOFDIR + '/evc' + aType + '.' + mstr + '.pkl'
        fid = open(fileevc,'r')
        evc[aType] = pkl.load(fid)
        fid.close()
        fileeva = EOFDIR + '/eva' + aType + '.' + mstr + '.pkl'
        fid = open(fileeva,'r')
        eva[aType] = pkl.load(fid)
        fid.close()
    
# For open sea
    print('open sea')
    for ii in range(nr0):
        ip = I0[ii]
        jp = J0[ii]
        inds = SUBtot[jp,ip]
        ireg = int(regs[jp,ip])-1
        isopen[ireg] = True
        pp = evc['open'][inds,:neof,:]
#        pp = EOF.variables['evc'].data[:neof,:,inds-1].copy()
        ap = eva['open'][inds,:neof]
#        ap = EOF.variables['eva'].data[:neof,inds-1].copy()
        var2D = varmap[jp,ip]
        if np.isnan(var2D) : print('Var NaN')
        lmbd = ap**2
        e2   = pp[:,0]**2
        sigma2 = np.sum(lmbd * e2)
        alpha = var2D/sigma2
        eofp[:,:ind200+1,ireg] = pp*(alpha**.25)
        eofa[:,ireg] = ap*(alpha**.25)

    print('coast ')
    pnoclass = 0
    for ii in range(nr2):
        ip = I2[ii]
        jp = J2[ii]
        inds = SUBtot[jp,ip]
        ireg = int(regs[jp,ip])-1
        indc = mapc[jp,ip]
        if indc==0:
            if isopen[ireg]==False:
                print('Point with not defined class')
                pnoclass += 1
#                import sys
#                sys.exit()
        else:
            cstr = '%.1i' %(indc)
            pp = evc['coast' + cstr][inds,:neof,:]
            ap = eva['coast' + cstr][inds,:neof]
            var2D = varmap[jp,ip]
            if np.isnan(var2D) : print('Var NaN')
            lmbd = ap[:neof]**2
            e2   = pp[:neof,0]**2
            sigma2 = np.sum(lmbd * e2)
            alpha = var2D/sigma2
            eofp[:,:ind200+1,ireg] = pp[:neof,:]*(alpha**.25)
            eofa[:,ireg] = ap[:neof]*(alpha**.25)

    nzero = np.sum(eofa[0,:]==0)
    if nzero<5:
        indz = np.where(eofa[0,:]==0)
        for iiz in indz:
            eofa[:,iiz] = eofa[:,iiz+1]
            eofp[:,:,iiz] = eofp[:,:,iiz+1]
    else: 
        print 'To much 0 in eof and eva EXIT'
        import sys
        sys.exit()
        

    fileout = OUTEOF + '/eof.' + mstr + '.nc'
    print('...Saving ' + fileout)
    write_eof(fileout,eofp,eofa)

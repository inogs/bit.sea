import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Eofs assignment for each grid point
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

    parser.add_argument(   '--grid3dvar', '-g',
                            type = str,
                            required = True,
                            help = '3dvar grid'
                            )

    parser.add_argument(   '--nleveof', '-n',
                            type = int,
                            required = True,
                            help = 'nleveof'
                            )

    return parser.parse_args()

args = argument()



import numpy as np
from bitsea.basins import V2
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
import scipy.io.netcdf as NC

from bitsea.commons.utils import addsep

INDIR = addsep(args.ineof)
OUTDIR = addsep(args.outdir)

BFMgrid = args.grid3dvar

meshmask = args.maskfile
TheMask = Mask(meshmask)

MM = NC.netcdf_file(BFMgrid,'r')
nregs = np.int(np.max(MM.variables['regs'].data))
regs = MM.variables['regs'].data.copy()
tgr = MM.variables['tmsk'].data[0,:,:].copy()
Lon = MM.variables['lon'].data.copy()
Lat = MM.variables['lat'].data.copy()

MM.close()


neof = 26
jpk,_,_ = TheMask.shape
nleveof = args.nleveof

# case = 'MONTHGRCLIM'

[JJ,II] = np.where(tgr==1)
ncells = II.shape[0]
reg_nomask = False
if nregs<ncells:
    # print ('Nmber of regions not coincident with tmask points')
    raise ValueError('Nmber of regions not coincident with tmask points')
if nregs>ncells:
    reg_nomask = True
    mm = (tgr==0) & (regs>1)
    indnomask = regs[mm]
    Nnomask = indnomask.shape[0]
    if Nnomask>3:
        raise ValueError('To much points with regs numebr but not in mask')
    IInomask = np.zeros(Nnomask,int)
    JJnomask = np.zeros(Nnomask,int)
    for iip in range(Nnomask):
        [JJnomask[iip],IInomask[iip]] = np.where(regs==int(indnomask[iip]))

LISTsub_all = [sub.name for sub in V2.Pred]

adrseen = False
LISTsubind = []
for isub,sub in enumerate(V2.Pred):
    if 'adr2' in sub.name:
        subi = isub-1
        adrseen = True
    else:
        if not(adrseen):
            subi = isub
        else:
            subi = isub-1
    LISTsubind.append(subi)

subindexes = np.zeros_like(tgr,dtype=int)
for isub,sub in enumerate(V2.Pred):
    print (sub.name)
    submask = SubMask(sub,maskobject=TheMask).mask[0,:,:]
    subindexes[submask] = LISTsubind[isub]
    
def write_eof(outfile,eofp,eofa):
    ncOUT   = NC.netcdf_file(outfile,"w")

    ncOUT.createDimension('nreg',nregs)
    ncOUT.createDimension('nlev',jpk)
    ncOUT.createDimension('neof',neof)

    ncvar = ncOUT.createVariable('eva','f',('neof','nreg'))
    ncvar[:] = eofa
    ncvar = ncOUT.createVariable('evc','f',('neof','nlev','nreg'))
    ncvar[:] = eofp

    ncOUT.close()
    

for im in range(12):
    eofp = np.zeros((neof,jpk,nregs))
    eofa = np.zeros((neof,nregs))
    mm = '%02d' %(im+1)
    print (mm)
    eoffile = INDIR + 'eof.15.' + mm + '.nc'
    EOF = NC.netcdf_file(eoffile,'r')
    evc = EOF.variables['evc'].data[:,:,:].copy()
    eva = EOF.variables['eva'].data[:,:].copy()
    EOF.close()
    for iip in range(ncells):
        ip = II[iip]
        jp = JJ[iip]
        indsub = subindexes[jp,ip]
        eofp[:,:nleveof,iip] = evc[:,:,indsub]
        eofa[:,iip] = eva[:,indsub]
    if reg_nomask:
        for iip in range(Nnomask):
            ip = IInomask[iip]
            jp = JJnomask[iip]
            indusb = subindexes[jp,ip]
            eofp[:,:nleveof,int(indnomask[iip])] = evc[:,:,indsub]
            eofa[:,int(indnomask[iip])] = eva[:,indsub]

    fileout = OUTDIR +  '/eof.' + np.str(nregs) + '.' + mm + '.nc' 
    print (fileout)
    write_eof(fileout,eofp,eofa)
    
    

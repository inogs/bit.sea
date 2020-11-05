import numpy as np
from basins import V2
from commons.mask import Mask
from commons.submask import SubMask
import scipy.io.netcdf as NC

INDIR = 'EOF15sub_600/'

BFMgrid = 'BFM_grid.nc'

meshmask = '/gpfs/scratch/userexternal/ateruzzi/MASKS16corrected/meshmask.nc'
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
nleveof = 40

case = 'MONTHGRCLIM'

[JJ,II] = np.where(tgr==1)
ncells = II.shape[0]
if not(nregs==ncells):
    print 'Nmber of regions not coincident with tmask points'
    import sys
    sys.exit(0)

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
    print sub.name
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
    print mm
    eoffile = INDIR + 'eof.15.' + mm + '.nc'
    EOF = NC.netcdf_file(eoffile,'r')
    evc = EOF.variables['evc'].data[:,:,:].copy()
    eva = EOF.variables['eva'].data[:,:].copy()
    EOF.close()
    for iip in range(nregs):
        ip = II[iip]
        jp = JJ[iip]
        indsub = subindexes[jp,ip]
        eofp[:,:nleveof,iip] = evc[:,:,indsub]
        eofa[:,iip] = eva[:,indsub]

    fileout = 'OUTEOF/' + case + '/eof.' + np.str(nregs) + '.' + mm + '.nc' 
    write_eof(fileout,eofp,eofa)
    
    

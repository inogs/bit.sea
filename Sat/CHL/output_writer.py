import netCDF4
import numpy as np
import dom_dec



maskfile="CHL_1km_meshmask.nc"


D=netCDF4.Dataset(maskfile,'r')
tmask_glo = np.array(D.variables['tmask']).astype(np.bool)
D.close()

jpjglo,jpiglo = tmask_glo.shape

nproc_i = 4
nproc_j = 4
PROCESSES = np.arange(nproc_i*nproc_j)

CLIM = np.zeros((365,jpjglo,jpiglo), dtype=[('NUMB',np.int32), ('MEAN',np.float32),('STD',np.float32)])

for ip in PROCESSES:
    (jproc, iproc) = divmod(ip,nproc_i)
    inputfile= "CHL_clim_%d_%d.npy" %(jproc,iproc)
    I_start,I_end, J_start, J_end = dom_dec.dom_dec(iproc, jproc, jpiglo, jpjglo, nproc_i, nproc_j)
    tmask = tmask_glo[J_start:J_end, I_start:I_end]
    print "process " , ip, inputfile, iproc, jproc, I_start,I_end, J_start, J_end
    Nwaterpoints = tmask.sum()
    if (Nwaterpoints ==0) : continue

    local_clim = np.load(inputfile)
    print local_clim['MEAN'].max()
    
    CLIM[:,J_start:J_end,I_start:I_end] = local_clim



ii=CLIM['NUMB']==0
CLIM['MEAN'][ii]=-999
CLIM['STD' ][ii]=-999

ncOUT = netCDF4.Dataset('Climatology_CHL.nc','w')
ncOUT.createDimension('lon',jpiglo)
ncOUT.createDimension('lat',jpjglo)
ncOUT.createDimension('time',365)
ncvar=ncOUT.createVariable('Mean','f',('time','lat','lon'))
setattr(ncvar,'missing_value',-999.0)
ncvar[:]=CLIM['MEAN']
ncvar=ncOUT.createVariable('Std','f',('time','lat','lon'))
setattr(ncvar,'missing_value',-999.0)
ncvar[:]=CLIM['STD']
ncvar=ncOUT.createVariable('Numb','i',('time','lat','lon'))
setattr(ncvar,'missing_value',0)
ncvar[:]=CLIM['NUMB']

ncOUT.close()

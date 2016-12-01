import numpy as np
import os,time
import scipy.io.netcdf as NC
from commons.dataextractor import DataExtractor
import netCDF4

class filename_manager():
    def __init__(self,filename):
        sep = "."
        nsep = os.path.basename(filename).count(".")
        if nsep ==3:
            prefix, datestr, varname,_ = os.path.basename(filename).rsplit(sep)
            self.varname = varname
        if nsep ==2:
            prefix, datestr, _ = os.path.basename(filename).rsplit(sep)
        self.prefix = prefix
        self.datestr= datestr

    def search_in_postproc_files(self, filename_for_timelist,var, AGGREGATE_AVEDIR):
        '''
        Arguments:
        filename_from_timelist (str) : POSTPROC/AVE_FREQ_1/ave.20130101-12:00:00.nc
        If var is an aggregate one, we search first of all in AGGREGATE AVEDIR
        RD, a read_descriptor object
        '''
        
        file_try1 = filename_for_timelist
        prefix, datestr, _ = os.path.basename(filename).rsplit(".")
        
        file_try2 = AGGREGATE_AVEDIR + prefix + "." + datestr + "."  + var + ".nc"

        if os.path.exists(file_try1):
            d=netCDF4.Dataset(file_try1,'r')
            if var in d.variables:
                d.close()
                return file_try1
        if os.path.exists(file_try2):
            d=netCDF4.Dataset(file_try2,'r')
            if var in d.variables:
                d.close()
                return file_try2
        raise ValueError("File not found")        



    def search_in_model_files(self, filename_for_timelist,var,INPUT_AVEDIR,AGGREGATE_AVEDIR):
        '''
        Arguments:
        filename_from_timelist (str) : MODEL/AVE_FREQ_1/ave.20130101-12:00:00.N1p.nc

        Returns  
        POSTPROC/output/AVE_FREQ_1/TMP/ave.20130101-12:00:00.P_l.nc
        '''
        
        prefix, datestr, _,_ = os.path.basename(filename_for_timelist).rsplit(".")
        file_try1 = INPUT_AVEDIR      + prefix + "." + datestr + "." + var + ".nc"
        file_try2 = AGGREGATE_AVEDIR  + prefix + "." + datestr + "." + var + ".nc"
        
        if os.path.exists(file_try1):
            d=netCDF4.Dataset(file_try1,'r')
            if var in d.variables:
                d.close()
                return file_try1
        if os.path.exists(file_try2):
            d=netCDF4.Dataset(file_try2,'r')
            if var in d.variables:
                d.close()
                return file_try2
        raise ValueError("File not found")


    def get_filename(self, filename_for_timelist,var,INPUT_AVEDIR,AGGREGATE_AVEDIR):
        '''
        Arguments:
        filename_from_timelist (str) : MODEL/AVE_FREQ_1/ave.20130101-12:00:00.N1p.nc
                                       POSTPROC/AVE_FREQ_1/ave.20130101-12:00:00.nc
        Returns
        POSTPROC/AVE_FREQ_1/TMP/ave.20130101-12:00:00.nc
        POSTPROC/AVE_FREQ_1/TMP/ave.20130101-12:00:00.P_l.nc
        MODEL/AVE_FREQ_1/ave.20130101-12:00:00.N3n.nc
        '''
        if is_a_big_avefile(filename_for_timelist):
            return self.search_in_postproc_files(filename_for_timelist, var,AGGREGATE_AVEDIR)
        else:
            return self.search_in_model_files(filename_for_timelist, var, INPUT_AVEDIR,AGGREGATE_AVEDIR)
        
            
        

def MoreRecent(file1,file2):
    A=time.gmtime(os.path.getmtime(file1))
    B=time.gmtime(os.path.getmtime(file2))
    return A>B

def is_a_big_avefile(filename):
    return os.path.basename(filename).count('.')==2

def getfileForRead(N1pfile, var):
    '''Can change var name in a file name
    ave.20130101-12:00:00.N1p.nc --> ave.20130101-12:00:00.N3n.nc '''
    if is_a_big_avefile(N1pfile):
        return N1pfile
    else:
        F = filename_manager(N1pfile)
        return os.path.dirname(N1pfile) + "/" + F.prefix + "." + F.datestr + "." + var + ".nc"


def WriteAggregateAvefiles(mask, N1pfile,OUTDIR,VarDescriptor):
    '''
    Writes one file for each aggregate var
      reads with NetCDF4
      writes in NetCDF3
    '''
            
    jpk,jpj,jpi = mask.shape
    F = filename_manager(N1pfile)
        
    for var in VarDescriptor.AGGREGATE_DICT.keys():
        outfile = OUTDIR + F.prefix + "." + F.datestr + "." + var + ".nc"
        ncOUT=NC.netcdf_file(outfile,"w")
        setattr(ncOUT,"Convenctions","COARDS")
        ncOUT.createDimension('time',  1)
        ncOUT.createDimension('depth', jpk)
        ncOUT.createDimension('lat' ,  jpj)
        ncOUT.createDimension('lon'  , jpi)
        ncvar=ncOUT.createVariable('lon','f',('lon',))
        setattr(ncvar,"units"       ,"degrees_east")
        ncvar[:]=mask.xlevels[0,:]
        ncvar=ncOUT.createVariable('lat','f',('lat',))
        setattr(ncvar,"units"       ,"degrees_north")
        ncvar[:]=mask.ylevels[:,0]
        ncvar=ncOUT.createVariable('depth','f',('depth',))
        setattr(ncvar,"units"       ,"meter")
        ncvar[:]=mask.zlevels
        setattr(ncOUT.variables['lon'],"long_name","Longitude")    
        setattr(ncOUT.variables['lat'],"long_name","Latitude")
        
        ncvar=ncOUT.createVariable(var,'f',('time','depth','lat','lon'))
        junk = np.zeros((jpk,jpj,jpi),np.float32)
        for lvar in VarDescriptor.AGGREGATE_DICT[var]:
            avefile=getfileForRead(N1pfile, lvar)
            DE = DataExtractor(mask,avefile,lvar)   
            junk +=DE.values
        junk[~mask.mask] = 1.e+20
        ncvar[:]=junk
        setattr(ncvar,"long_name",var)
        setattr(ncvar,"missing_value",1.e+20) 
        del DE
        ncOUT.close()      
    

def WriteBigAve(Mask,N1pfile, outfile, VARS):
      
    if len(VARS)==0:
        print "No variables in archive list"
        return
            
    nc=NC.netcdf_file(N1pfile,"r");
    DIMS=nc.dimensions;


    
    ncOUT=NC.netcdf_file(outfile,"w")
    setattr(ncOUT,"Convenctions","COARDS")
    setattr(ncOUT,"DateStart",nc.DateStart)
    setattr(ncOUT,"Date__End",nc.Date__End)
        
    
    for dimName,dimValue in DIMS.items():
        ncOUT.createDimension(dimName,dimValue)
    
    for var in ['lon','lat','depth']:
        ncvar=ncOUT.createVariable(var,'f',(var,))
        ncvar[:]=nc.variables[var].data.copy()
        setattr(ncvar,"actual_range",nc.variables[var].actual_range)
        setattr(ncvar,"units"       ,nc.variables[var].units)
    nc.close()        
    
    setattr(ncOUT.variables['lon'],"long_name","Longitude")    
    setattr(ncOUT.variables['lat'],"long_name","Latitude")

    
    for var in VARS:
        avefile=getfileForRead(N1pfile, outfile, var)        
        ncIN = NC.netcdf_file(avefile,"r")    
        ncvar=ncOUT.createVariable(var,'f',('time','depth','lat','lon'))
        ncvar[:]=ncIN.variables[var].data.copy()
        setattr(ncvar,"long_name",var)
        setattr(ncvar,"missing_value",1.e+20)        
        ncIN.close()

    ncOUT.close()    



def writeChlSup(avefile, chlfile, chlvar):
    '''
    chlvar can be P_l or P_i
    '''
    AVE=NC.netcdf_file(avefile,"r")
    CHL=NC.netcdf_file(chlfile,"w")
    fillValue=1.e+20
    
    a=AVE.variables[chlvar].data[0,0,:,:].copy()      
    a[a>fillValue]=fillValue
    
    DIMS=AVE.dimensions
    CHL.createDimension('time',None)
    CHL.createDimension('lon',DIMS['lon'])
    CHL.createDimension('lat',DIMS['lat'])
    ncvar=CHL.createVariable("lchlm",'f',('lat','lon'))
    ncvar[:]=a
    setattr(ncvar,"missing_value",fillValue)
    AVE.close()
    CHL.close()

if __name__ == "__main__" :
    import read_descriptor

    
    RD=read_descriptor.read_descriptor('VarDescriptor_2.xml')
        
    
    var='P_i'
    INPUT_AVEDIR = "MODEL/AVE_FREQ_1/"
    AGGREGATE_AVEDIR="POSTPROC/TMP/"
    
    avelist=["ave.20130101-12:00:00.N1p.nc","ave.20130102-12:00:00.N1p.nc","ave.20130103-12:00:00.N1p.nc"]
    avelist=["ave.20130101-12:00:00.nc","ave.20130102-12:00:00.nc","ave.20130103-12:00:00.nc"]
    avelist=[INPUT_AVEDIR + f for f in avelist]
    
    filename=avelist[0]
    F=filename_manager(filename)
    print F.get_filename(filename, var,RD, INPUT_AVEDIR,AGGREGATE_AVEDIR)
    

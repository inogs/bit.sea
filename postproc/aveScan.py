import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    
    With -s flag, creates statistics on ave files contained in tmpdir, i.e. already aggregated.
    Statistics are defined on  on subbasins, then a submask is needed (as environment variable).
    Statistical profiles are [mean, std, p25, p50, p75], of each level, weighted on cell area over a subbasin.
    If a variable is defined in Some_Statistics group, percentiles are skipped
    Integrals are [mean, std, p25, p50, p75], for each depth in DEPTHlist, weighted on cell volume over a subbasin.
    DEPTHlist is provided and defined in maskload (shallow = [0,200m], deep=[200,bottom]) 
    For variables in Some_Statistics group, percentiles are skipped and Volume integrals 
    are based on statistical profiles.
    For the only ppn var (if exists) a different kind of column integral is calculaded. 
    
    After that, there are point profiles, at locations defined in --pointlist argument.
    Without this argument, extraction of point profiles is skipped.
    ''')
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = 'The directory which will contain INTEGRALS/, PROFILES/ and STAT_PROFILES/')

    parser.add_argument(   '--tmpdir', '-t',
                                type = str,
                                required =True,
                                help = '''
                                Path of the temporary files tmp*.profiles. or tmp*integrals.nc
                                ''')
        
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required =True,
                                help = '''
                                Path of the 1st input directory. Usually it is the MODEL/AVE_FREQ_1/ directory produced by ogstm.xx.
                                Here variable are considered as native.
                                ''')
    
    parser.add_argument(   '--aggregatedir', '-a',
                                type = str,
                                required =True,
                                help = '''
                                Path of the 2nd input directory. Usually it is the TMP directory produced by postproc.py.
                                Here variables are considered as already aggregated.
                                ''')

    
    parser.add_argument(   '--avelist',"-l",
                                type = str,
                                required = True,
                                help = 'ave*nc')
    parser.add_argument(   '--var',"-f",
                                type = str,
                                required = False,
                                help = ''' 
                                var name, like P_l useful to filter the file list. 
                                No filtering is done if it is not provided.
                                 ''')
    
    parser.add_argument(   '--descriptor',"-d",
                                type = str,
                                default = "VarDescriptor_1.xml",
                                help = 'VarDescriptor_1.xml, or the complete path') 
    
    parser.add_argument(   '-s',action='store_true',
                                help = '''Activates statistics calculation
                                ''') 
    parser.add_argument(   '--pointlist',"-p",
                                type = str,
                                help = '''Path of the text file listing the the points where extract point profiles''')             
    parser.add_argument(   '--ionames',
                                type = str,
                                default="IOnames.xml",
                                help = '''Path of IOnames file''')

    return parser.parse_args()




args = argument()

import scipy.io.netcdf as NC
import netCDF4 as NC4
from commons.dataextractor import DataExtractor
import glob
import os
import numpy as np
import read_descriptor
import IOnames as IOname
from maskload import *
import GB_lib
from commons.utils import addsep


INPUT_AVEDIR = addsep(args.inputdir)
AGGREGATE_AVEDIR = addsep(args.aggregatedir)
TMPDIR       = addsep(args.tmpdir)
BASEDIR      = addsep(args.outdir)
ionamesfile  = args.ionames
IOnames      = IOname.IOnames(ionamesfile)
filtervar    = args.var

rank   = 0
nranks = 1
isParallel = False
    
if rank==0 : print("list      " ,args.avelist, "filtered by", filtervar)
def wp(data, wt, percentiles):
    """Compute weighted percentiles. 
    If the weights are equal, this is the same as normal percentiles. 
    Elements of the C{data} and C{wt} arrays correspond to 
    each other and must have equal length (unless C{wt} is C{None}). 
    @param data: The data. 
    @type data: A L{np.ndarray} array or a C{list} of numbers. 
    @param wt: How important is a given piece of data. 
    @type wt: C{None} or a L{np.ndarray} array or a C{list} of numbers. 
    All the weights must be non-negative and the sum must be 
    greater than zero. 
    @param percentiles: what percentiles to use. (Not really percentiles, 
    as the range is 0-1 rather than 0-100.) 
    @type percentiles: a C{list} of numbers between 0 and 1. 
    @rtype: [ C{float}, ... ] 
    @return: the weighted percentiles of the data. 
    """ 
    assert np.greater_equal(percentiles, 0.0).all(), "Percentiles less than zero" 
    assert np.less_equal(percentiles, 1.0).all(), "Percentiles greater than one"
    data = np.asarray(data) 
    assert len(data.shape) == 1 
    if wt is None:
        wt = np.ones(data.shape, np.float)     
    else: 
        wt = np.asarray(wt, np.float) 
        assert wt.shape == data.shape 
        assert np.greater_equal(wt, 0.0).all(), "Not all weights are non-negative." 
    assert len(wt.shape) == 1
    n = data.shape[0] 
    assert n > 0 
    i = np.argsort(data) 
    sd = np.take(data, i, axis=0) 
    sw = np.take(wt, i, axis=0) 
    aw = np.add.accumulate(sw) 
    if not aw[-1] > 0: 
        raise ValueError("Nonpositive weight sum")
    w = (aw-0.5*sw)/aw[-1] 
    spots = np.searchsorted(w, percentiles) 
    o = [] 
    for (s, p) in zip(spots, percentiles): 
        if s == 0: 
            o.append(sd[0]) 
        elif s == n:
            o.append(sd[n-1]) 
        else: 
            f1 = (w[s] - p)/(w[s] - w[s-1]) 
            f2 = (p - w[s-1])/(w[s] - w[s-1]) 
            assert f1>=0 and f2>=0 and f1<=1 and f2<=1 
            assert abs(f1+f2-1.0) < 1e-6 
            o.append(sd[s-1]*f1 + sd[s]*f2) 
    return o 

#nstat = 5 ; StatDescr   = "Mean, Std, p25, p50, p75"
nstat = 9  ; StatDescr   = "Mean, Std, min, p05, p25, p50, p75, p95, max"
def CoreStatistics(Conc, Weight):
    Statistics      = np.zeros((nstat),np.float32)
    
    nans = Conc == 1.e+20
    Conc = Conc[~nans]
    Weight = Weight[~nans]
    if Conc.nbytes ==0 :
            return Statistics   
    
    Weight_sum      = Weight.sum()
    Mass            = (Conc * Weight).sum()
    Weighted_Mean   = Mass/Weight_sum        
    Weighted_Std    = ((Conc - Weighted_Mean)**2*Weight).sum()/Weight_sum
    
    Statistics[0]   = Weighted_Mean
    Statistics[1]   = np.sqrt(Weighted_Std)
    #Statistics[2:5] = wp(Conc,Weight,[.25, .50,.75]) #for nstat=5

    Statistics[2]   = Conc.min()
    Statistics[3:8] = wp(Conc,Weight,[.05, .25, .50,.75, .95])
    Statistics[8]   = Conc.max()
    
    return Statistics



def CoreStatistics_noSort(Conc,Weight):    
    Statistics      = np.zeros((nstat),np.float32)
    if Conc.nbytes ==0 :
            return Statistics   
    
    Weight_sum      = Weight.sum()
    if  Weight_sum == 0:
        return Statistics
    Mass            = (Conc * Weight).sum()
    Weighted_Mean   = Mass/Weight_sum
    Weighted_Std    = ((Conc - Weighted_Mean)**2*Weight).sum()/Weight_sum      
    Statistics[0]   = Weighted_Mean
    Statistics[1]   = np.sqrt(Weighted_Std)
    Statistics[2]   = Conc.min()
    Statistics[8]   = Conc.max()
    return Statistics

def Population_Statistics(averages, stds, Weight):
    Statistics      = np.zeros((2),np.float32)
    if averages.nbytes ==0 :
            return Statistics   
    
    Weight_sum      = Weight.sum()
    if  Weight_sum == 0:
        return Statistics
        
    Mass       = (averages * Weight).sum()
    Pop_Mean   = Mass/Weight_sum
    #Pop_Std    = np.sqrt(  ( (averages**2 + stds)*Weight).sum()/Weight_sum - Pop_Mean**2)
    Pop_Std    = (Weight/Weight_sum * stds).sum()
    Statistics[0]   = Pop_Mean
    Statistics[1]   = Pop_Std
    return Statistics


def PointProfiles(varname):
    nCruise = len(MeasPoints)
    Profiles = np.zeros((nCruise,jpk),np.float32)
    
    for c in range(nCruise):
        point=MeasPoints[c]
        name = point['Name']
        j    = point['j']
        i    = point['i']        
        Profiles[c,:] = VAR[:,j,i]
        
    return Profiles
    


def ProfilesStats(Mask, flagSort):
# VAR and area are read from main        
    Statistics   = np.zeros((jpk,nstat),np.float32)    
    
    if flagSort:    
        for k in range(jpk):
            V =  VAR[k,:,:]
            M = Mask[k,:,:]     
            Statistics[k,:] = CoreStatistics(V[M],area[M])
    else:
        for k in range(jpk):
            V =  VAR[k,:,:]
            M = Mask[k,:,:]                    
            Statistics[k,:] = CoreStatistics_noSort(V[M],area[M])                
    
    return Statistics
    
    


     

def Vol_Integrals(Mask, flagSort=True):

    Statistics   = np.zeros((ndepth,nstat),np.float32)    
     
    if flagSort:
        for idepth in range(len(DEPTHlist)):
            depth = DEPTHlist[idepth]
            mask  = Mask & DEPTH[depth]
            Statistics[idepth,:] = CoreStatistics(VAR[mask],Volume[mask])
    else:
        for idepth in range(len(DEPTHlist)):
            depth = DEPTHlist[idepth]
            mask  = Mask & DEPTH[depth]
            Statistics[idepth,:] = CoreStatistics_noSort(VAR[mask],Volume[mask])

    return Statistics

def Vol_Integrals_On_Profile(Mask, profileM, profileSTD):

    Statistics   = np.zeros((ndepth,2),np.float32)    
    Level_Volume = np.zeros((jpk),np.float32)
     
    for idepth in range(len(DEPTHlist)):
        depth = DEPTHlist[idepth]
        mask  = Mask & DEPTH[depth]        
        for k in range(jpk):
            Vol=Volume[k,:,:]
            sm =  mask[k,:,:]
            Level_Volume[k] = Vol[sm].sum()
        Statistics[idepth,:]= Population_Statistics(profileM, profileSTD, Level_Volume)         

    return Statistics
    

def getColumnIntegrals(var, objfileI):
    INTEGRALS    = np.zeros((nsub,ncoast,ndepth,nstat),np.float32)
    iDepth = 0
    iStat  = 0
    for isub in range(len(SUBlist)):
        sub   = SUBlist[isub]
        m = SUB(sub)
        for icoast in range(len(COASTNESS_LIST)):
            
            coast = COASTNESS_LIST[icoast]
            smask = m & COASTNESS[coast]
            Cdz   = VAR*dZ*smask
            Int_Cdz = np.nansum(Cdz,axis=0)
            if np.any(smask[0,:,:]):
                INTEGRALS[isub,icoast,iDepth,iStat] = Int_Cdz[smask[0,:,:]].mean(axis=None)
            else:
                print("getColumnIntegrals: no data found in ", sub, coast)
    ncvar    = objfileI.createVariable(var, 'f', ('sub','coast','depth','stat'))
    ncvar[:] = INTEGRALS  
      
    
def getSomeStatistics(var,objfileP,objfileI):
    #print datetime.datetime.now()
    PROFILES     = np.zeros((nsub,ncoast,jpk   ,nstat),np.float32)
    INTEGRALS    = np.zeros((nsub,ncoast,ndepth,nstat),np.float32)
    
    
    for isub in range(len(SUBlist)):
        sub   = SUBlist[isub]
        m = SUB(sub)
        
        for icoast in range(len(COASTNESS_LIST)):
            
            coast = COASTNESS_LIST[icoast]
            smask = m & COASTNESS[coast]
            PROFILES [isub,icoast,:,:] = ProfilesStats(smask, flagSort=False)    
            #INTEGRALS[isub,icoast,:,:]= Vol_Integrals_On_Profile(smask, PROFILES [isub,icoast,:,0])
        smask= m & COASTNESS['open_sea']
        profile_mean = PROFILES[isub,1,:,0]
        profile_std  = PROFILES[isub,1,:,1]
        INTEGRALS[isub,1,:,:2]= Vol_Integrals_On_Profile(smask, profile_mean, profile_std)
        
    ncvar    = objfileP.createVariable(var  ,'f',('sub','coast','z'    ,'stat'))    
    ncvar[:] = PROFILES        
    ncvar    = objfileI.createVariable(var, 'f', ('sub','coast','depth','stat'))
    ncvar[:] = INTEGRALS        
        
 
    
    
def getAllStatistics(var,objfileP,objfileI):
    #here VAR can be read from main   
    #print datetime.datetime.now()
    
    PROFILES = np.zeros((nsub,ncoast,jpk   ,nstat),np.float32)
    INTEGRALS =np.zeros((nsub,ncoast,ndepth,nstat),np.float32)
    
    for isub in range(len(SUBlist)):
        sub   = SUBlist[isub]
        m = SUB(sub)
        
        for icoast in range(len(COASTNESS_LIST)):
            
            coast = COASTNESS_LIST[icoast]
            smask = m & COASTNESS[coast]
            PROFILES [isub,icoast,:,:] = ProfilesStats(smask, flagSort=True)
            INTEGRALS[isub,icoast,:,:] = Vol_Integrals(smask, flagSort=True)
        
        
    
    ncvar    = objfileP.createVariable(var  ,'f',('sub','coast','z'    ,'stat'))    
    ncvar[:] = PROFILES  
    ncvar    = objfileI.createVariable(var, 'f', ('sub','coast','depth','stat'))
    ncvar[:] = INTEGRALS

def ProfilesStats2D(Mask, flagSort):
# VAR and area are read from main        
    Statistics   = np.zeros((jpk,nstat),np.float32)    
    
    if flagSort:    
        for k in range(1):
            V =  VAR[k,:,:]
            M = Mask[k,:,:]
            
            Statistics[k,:] = CoreStatistics(V[M],area[M])
    else:
        for k in range(1):
            V =  VAR
            M = Mask[k,:,:]                    
            Statistics[k,:] = CoreStatistics_noSort(V[M],area[M])                
    
    return Statistics
    

def getAllStatistics2D(var,objfileP):
    #here VAR can be read from main   
    
    PROFILES = np.zeros((nsub,ncoast,jpk   ,nstat),np.float32)
    
    for isub in range(len(SUBlist)):
        sub   = SUBlist[isub]
        m = SUB(sub)
        
        for icoast in range(len(COASTNESS_LIST)):
            
            coast = COASTNESS_LIST[icoast]
            smask = m & COASTNESS[coast]
            PROFILES [isub,icoast,:,:] = ProfilesStats2D(smask, flagSort=True)                    
        
    
    ncvar    = objfileP.createVariable(var  ,'f',('sub','coast','z'    ,'stat'))    
    ncvar[:] = PROFILES  

    
def dumpPointProfiles(var):
    
# now VAR has to be modified
    tmp_Pprofiles = TMPDIR + "tmp." + datestr + "." + var + ".profiles.nc"
    ncOUT_Pprofiles = NC.netcdf_file(tmp_Pprofiles,"w")
    ncOUT_Pprofiles.createDimension("Ncruise"   ,nCruise)
    ncOUT_Pprofiles.createDimension("z"         ,jpk)            
    global VAR
    #VAR = VAR.copy()
    VAR[~tmask]=np.NaN    
    P=PointProfiles(var)
    ncvar    = ncOUT_Pprofiles.createVariable(var  ,'f',('Ncruise','z'))    
    ncvar[:] = P
    ncOUT_Pprofiles.close()       

def create_tmp_headers(datestr,var):
    tmp__profiles = TMPDIR + "tmp." + datestr + "." + var + ".stat_profiles.nc"
    tmp_integrals = TMPDIR + "tmp." + datestr + "." + var + ".vol_integrals.nc"     
       
    ncOUT__profiles = NC.netcdf_file(tmp__profiles,"w")
    ncOUT_integrals = NC.netcdf_file(tmp_integrals,"w")    
    
    ncOUT__profiles.createDimension("sub"       ,len(SUBlist))
    ncOUT__profiles.createDimension("coast"     ,ncoast)
    ncOUT__profiles.createDimension("z"         ,jpk )
    ncOUT__profiles.createDimension("stat"      ,nstat )
    
    ncOUT_integrals.createDimension("sub"       ,len(SUBlist))
    ncOUT_integrals.createDimension("coast"     ,ncoast)
    ncOUT_integrals.createDimension("depth"     ,ndepth)
    ncOUT_integrals.createDimension("stat"      ,nstat )
    
    return ncOUT__profiles, ncOUT_integrals      

def create_ave_headers(datestr):
    '''
    Generates output files, in INTEGRALS and STAT_PROFILES directories
    Files refer to a time frame and will contain all the variables
    Works in append mode, if the file name exists.
    '''
    
    ave__profiles = OUTPUT_DIR_MPR + IOnames.Output.prefix + datestr + IOnames.Output.suffix + ".stat_profiles.nc"
    ave_integrals = OUTPUT_DIR_INT + IOnames.Output.prefix + datestr + IOnames.Output.suffix + ".vol_integrals.nc" 

    if os.path.exists(ave__profiles):
        ncOUT__profiles=NC4.Dataset(ave__profiles,'a')
        print("appending in ", ave__profiles)
    else:
        ncOUT__profiles = NC4.Dataset(ave__profiles,'w')
        ncOUT__profiles.createDimension("sub"       ,len(SUBlist))
        ncOUT__profiles.createDimension("coast"     ,ncoast)
        ncOUT__profiles.createDimension("z"         ,jpk )
        ncOUT__profiles.createDimension("stat"      ,nstat )        
        setattr(ncOUT__profiles,"sub___list"  ,  SubDescr[:-2])
        setattr(ncOUT__profiles,"coast_list"  ,CoastDescr[:-2])
        setattr(ncOUT__profiles,"stat__list"  , StatDescr)

    if os.path.exists(ave_integrals):
        ncOUT_integrals=NC4.Dataset(ave_integrals,'a')
        print("appending in ", ave_integrals)
    else:
        ncOUT_integrals=NC4.Dataset(ave_integrals,'w')
        ncOUT_integrals.createDimension("sub"       ,len(SUBlist))
        ncOUT_integrals.createDimension("coast"     ,ncoast)
        ncOUT_integrals.createDimension("depth"     ,ndepth)
        ncOUT_integrals.createDimension("stat"      ,nstat )
    
        setattr(ncOUT_integrals,"sub___list"  ,  SubDescr[:-2])    
        setattr(ncOUT_integrals,"coast_list"  ,CoastDescr[:-2])
        setattr(ncOUT_integrals,"depth_list"  ,DepthDescr[:-2])
        setattr(ncOUT_integrals,"stat__list"  , StatDescr)
    
    return ncOUT__profiles,  ncOUT_integrals

def create_ave_pp_header(datestr):
    '''
    Generates output files, in PROFILES directory
    Files refer to a time frame and will contain all the variables
    Works in append mode, if the file name exists.'''
    
    ave_Pprofiles = OUTPUT_DIR_PRO + IOnames.Output.prefix + datestr + ".profiles.nc"
    if os.path.exists(ave_Pprofiles):
        ncOUT_Pprofiles =  NC.netcdf_file(ave_Pprofiles,"a")
    else:
        ncOUT_Pprofiles = NC.netcdf_file(ave_Pprofiles,"w")
        ncOUT_Pprofiles.createDimension("Ncruise"   ,nCruise)
        ncOUT_Pprofiles.createDimension("z"         ,jpk)
        setattr(ncOUT_Pprofiles,"CruiseIndex",CruiseDescr)
    return ncOUT_Pprofiles

doPointProfiles = False; 
if args.pointlist: 
    doPointProfiles=True;
    MeasPoints = read_Positions_for_Pointprofiles(args.pointlist) 
    nCruise         = len(MeasPoints)
    CruiseDescr =""
    for i in MeasPoints['Name']: CruiseDescr+=str(i) + ", "

doStatistics    = args.s
#nstat           =  5
ncoast          = len(COASTNESS_LIST)
ndepth          = len(DEPTHlist)
nTrans          = 51
nsub            = len(SUBlist)



SubDescr    ="" 
CoastDescr  =""
DepthDescr  =""



for i in SUBlist           : SubDescr   +=str(i) + ", "
for i in COASTNESS_LIST    : CoastDescr +=str(i) + ", "
for i in DEPTHlist         : DepthDescr +=str(i) + ", "


RD = read_descriptor.read_descriptor(args.descriptor)
NATIVE_VARS    = RD.NATIVE_VARS
AGGREGATE_VARS = RD.AGGREGATE_VARS
SOME_VARS      = RD.SOME_VARS

OUTPUT_DIR_INT  = BASEDIR + 'INTEGRALS/'
OUTPUT_DIR_MPR  = BASEDIR + 'STAT_PROFILES/'
OUTPUT_DIR_PRO  = BASEDIR + 'PROFILES/'
OUTPUT_DIR_PPN  = BASEDIR + 'INTEGRALS/PPN/'

if rank==0 : 
    RD.printStatus() 
    print("INPUT_AVEDIR", INPUT_AVEDIR)
 
    for DIR in [OUTPUT_DIR_INT, OUTPUT_DIR_MPR, OUTPUT_DIR_PPN, TMPDIR]:
        if doStatistics:    os.system("mkdir -p " + DIR)
    for DIR in [OUTPUT_DIR_PRO, TMPDIR]:
        if doPointProfiles: os.system("mkdir -p " + DIR)
    print(INPUT_AVEDIR + args.avelist)
    
aveLIST = glob.glob(INPUT_AVEDIR + args.avelist)
if not filtervar is None:
    aveLIST=[f for f in aveLIST if filtervar in os.path.basename(f) ]
aveLIST.sort()

VARLIST=[]
VARLIST.extend(NATIVE_VARS)
VARLIST.extend(AGGREGATE_VARS)
VARLIST.extend(SOME_VARS)
var_dim =['3D']*len(VARLIST)
VARLIST.extend(RD.vars2D)
var_dim.extend(['2D']*len(RD.vars2D))

nvars  = len(VARLIST)
nFiles = len(aveLIST)

assert nFiles > 0

PROCESSES = np.arange(nFiles * nvars)


if isParallel : comm.Barrier() 


for ip in PROCESSES[rank::nranks]:
    (ifile, ivar) = divmod(ip,nvars)
    var     = VARLIST[ivar]
    avefile =aveLIST[ifile]
    dim     = var_dim[ivar]

    F = GB_lib.filename_manager(avefile)
    datestr=F.datestr# os.path.basename(avefile)[IOnames.Input.date_startpos:IOnames.Input.date_endpos]
    if doStatistics:
        ncOUT__profiles,ncOUT_integrals = create_tmp_headers(datestr,var)

    filename = F.get_filename(avefile, var,INPUT_AVEDIR,AGGREGATE_AVEDIR)
    print("rank %03d scans %s on %s" %(rank,var,os.path.basename(filename)))
    vartoread = var

    if ('RST' in os.path.basename(filename)) and (not('before' in os.path.basename(filename))):
        vartoread = 'TRN' + var

    if var_dim [ivar] == '3D':
        VAR  = DataExtractor(TheMask,filename,vartoread,dimvar=3).values
    else:
        VAR        = np.zeros((jpk,jpj,jpi),np.float32)
        VAR[0,:,:] = DataExtractor(TheMask,filename,vartoread,dimvar=2).values
    
    if doStatistics:
        if var_dim [ivar] == '3D':
            if var in SOME_VARS:
                getSomeStatistics(var, ncOUT__profiles, ncOUT_integrals)
            else:
                getAllStatistics(var, ncOUT__profiles, ncOUT_integrals)
        if var_dim[ivar] == '2D' :
            getAllStatistics2D(var, ncOUT__profiles)
    if doPointProfiles : dumpPointProfiles(var)

    if doStatistics:
        ncOUT__profiles.close()
        ncOUT_integrals.close()
        
    
if isParallel : comm.Barrier()    
if rank == 0: print("RICOSTRUZIONE FILES")

for avefile in aveLIST[rank::nranks]:
    F = GB_lib.filename_manager(avefile)
    datestr=F.datestr #os.path.basename(avefile)[IOnames.Input.date_startpos:IOnames.Input.date_endpos]

    if doPointProfiles: ncOUT_Pprofiles = create_ave_pp_header(datestr)
    if doStatistics:    ncOUT__profiles,ncOUT_integrals = create_ave_headers(datestr)

    
    for ivar, var in enumerate(VARLIST):
        tmp__profiles = TMPDIR + "tmp." + datestr + "." + var + ".stat_profiles.nc"
        tmp_integrals = TMPDIR + "tmp." + datestr + "." + var + ".vol_integrals.nc"
        tmp_Pprofiles = TMPDIR + "tmp." + datestr + "." + var + ".profiles.nc"
        
        if doStatistics: 
            ncIN__profiles = NC.netcdf_file(tmp__profiles,"r")
            ncIN_integrals = NC.netcdf_file(tmp_integrals,"r")
            if var in ncOUT__profiles.variables.keys(): #existing variable
                ncvar=ncOUT__profiles.variables[var]
            else:
                ncvar    = ncOUT__profiles.createVariable(var  ,'f',('sub','coast','z'    ,'stat'))

            ncvar[:] = ncIN__profiles.variables[var].data.copy()

            if var_dim[ivar] == '3D':
                if var in ncOUT_integrals.variables.keys(): #existing variable
                    ncvar = ncOUT_integrals.variables[var]
                else:
                    ncvar    = ncOUT_integrals.createVariable(var, 'f', ('sub','coast','depth','stat'))
                ncvar[:] = ncIN_integrals.variables[var].data.copy()
            ncIN__profiles.close()
            ncIN_integrals.close()
        
        if doPointProfiles:
            ncIN_Pprofiles = NC.netcdf_file(tmp_Pprofiles,"r")
            ncvar    = ncOUT_Pprofiles.createVariable(var  ,'f',('Ncruise','z'))
            ncvar[:] = ncIN_Pprofiles.variables[var].data.copy()       
            ncIN_Pprofiles.close()
                

    if doPointProfiles: ncOUT_Pprofiles.close()
    if doStatistics: 
        ncOUT__profiles.close()
        ncOUT_integrals.close()
    
    
   
if rank==0 : print("RICOSTRUZIONE FINITA")
if isParallel: comm.Barrier() 

var = 'ppn'
if var not in VARLIST or not doStatistics : sys.exit()

for avefile in aveLIST[rank::nranks]:
    print(avefile)
    F=GB_lib.filename_manager(avefile)
    filename = F.get_filename(avefile, var,INPUT_AVEDIR,AGGREGATE_AVEDIR)
    VAR  = DataExtractor(TheMask,filename,var,dimvar=3).values 
    
    ave_integrals = OUTPUT_DIR_PPN + "ave." + F.datestr + ".col_integrals.nc"
    ncOUT_integrals = NC.netcdf_file(ave_integrals,"w")
    ncOUT_integrals.createDimension("sub"       ,len(SUBlist))
    ncOUT_integrals.createDimension("coast"     ,ncoast)
    ncOUT_integrals.createDimension("depth"     ,ndepth)
    ncOUT_integrals.createDimension("stat"      ,nstat )
    
    getColumnIntegrals(var, ncOUT_integrals)
    
    setattr(ncOUT_integrals,"sub___list"  ,  SubDescr[:-2])    
    setattr(ncOUT_integrals,"coast_list"  ,CoastDescr[:-2])
    setattr(ncOUT_integrals,"depth_list"  ,DepthDescr[:-2])
    setattr(ncOUT_integrals,"stat__list"  , StatDescr    )
    ncOUT_integrals.close()


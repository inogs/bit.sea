import numpy as np
import netCDF4
import os
from commons.utils import addsep
from Float.oxygen_saturation import oxy_sat

class checkreport():
    def __init__(self, linestr, nObs,nExcl,reason, depthexc):
        self.line = linestr
        self.reason=np.nan
        self.nObservations = nObs
        self.nExcluded =nExcl
        self.reason = reason
        self.depthexc = depthexc
    def __repr__(self):
        return "CheckReport object: %d Observations %d excluded " %(self.nObservations, self.nExcluded)


class check():
    def __init__(self,OUTDIR, verboselevel=1):
        '''
        if verboselevel = 1, a NetCDF will be printed out
        else, no files will be dumped
        '''
        if OUTDIR=="":
            self.outdir=""
        else:
            self.outdir=addsep(OUTDIR)
            os.system("mkdir -p " + OUTDIR)
        self.verboselevel=verboselevel

    def nitrate_check(self, model, ref, depth, p):
    
        bad=(np.isnan(model)) | (np.isnan(ref))
        MODEL = model[~bad]
        REF   = ref [~bad]
        DEPTH= depth[~bad]
        nP = len(DEPTH)
    
        mydiff = np.abs(MODEL-REF)
        flag1_array = (mydiff > 1)  & (DEPTH<=50)
        flag2_array = (mydiff > 2)  & (DEPTH>=250) & (DEPTH<=600)
        flag1 = flag1_array.sum() > 5
        flag2 = flag2_array.sum() > 5
    
        flag3_array = (REF > 3)  & (DEPTH<=15)
        flag3 = flag3_array.sum() > 0
    
    
        nobs = DEPTH.shape[0]
        nexcl = 0
        firstdepth = np.nan
        if flag2:
            depthind = np.where(flag2_array)[0][0]
            firstdepth = DEPTH[depthind]
            nexcl = np.sum(DEPTH>=firstdepth)
        elif (flag1 | flag3):
            nexcl = nobs
            
    
        line=""
        flag=np.nan
        if ( flag1  |  flag2 | flag3 ):
            if flag1: flag=1
            if flag2: flag=2
            if flag3: flag=-1
            if (flag1  & flag2) : flag=3
            line="%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %( p._my_float.wmo, p.time.strftime("%Y%m%d"), p.lon, p.lat , nobs, nexcl, flag )
        
            FLAG = np.zeros((nP), np.int)
            FLAG[flag1_array] =1
            FLAG[flag2_array] =2
            FLAG[flag3_array] =-1
        
            outncfile="%s%s_%s.nc"  %(self.outdir + 'N3n.', p.time.strftime("%Y%m%d"), p._my_float.wmo )
            if self.verboselevel==1:
                print "check dumps "  + outncfile
                ncOUT = netCDF4.Dataset(outncfile,'w')
                ncOUT.createDimension('depth',nP)
                ncvar = ncOUT.createVariable('model','f',('depth', ))
                ncvar[:]=MODEL
                ncvar = ncOUT.createVariable('float','f',('depth', ))
                ncvar[:]=REF
                ncvar = ncOUT.createVariable('flag' ,'i',('depth', ))
                ncvar[:]=FLAG
                setattr(ncOUT, 'longitude', p.lon)
                setattr(ncOUT, 'latitude' , p.lat)
                ncOUT.close()
        CR = checkreport(line, nobs,nexcl, flag, firstdepth)
        return CR

    def chlorophyll_check(self, model, ref, depth, p):
        
        bad=(np.isnan(model)) | (np.isnan(ref))
        MODEL = model[~bad]
        REF   = ref [~bad]
        DEPTH= depth[~bad]
        nP = len(DEPTH)
        
        mydiff = np.abs(MODEL-REF)
        flag1_array = (mydiff > 2)  & (DEPTH<=200)
        flag1 = flag1_array.sum() > 5
        # flag2 = flag2_array.sum() > 5
        nexcl = 0
    
        line=""
        flag=np.nan
        if ( flag1) :#  |  flag2): 
            if flag1: flag=1
            # if flag2: flag=2
            # if (flag1  & flag2) : flag=3
            line="%s\t%s\t%s\t%s\t%s\n" %( p._my_float.wmo, p.time.strftime("%Y%m%d"), p.lon, p.lat , flag)
            nexcl=nP
            FLAG = np.zeros((nP), np.int)
            FLAG[flag1_array] =1
            # FLAG[flag2_array] =2
            if self.verboselevel ==1 :
                outncfile="%s%s_%s.nc"  %(self.outdir + 'P_l.', p.time.strftime("%Y%m%d"), p._my_float.wmo )
                ncOUT = netCDF4.Dataset(outncfile,'w')
                ncOUT.createDimension('depth',nP)

                ncvar = ncOUT.createVariable('model','f',('depth', ))
                ncvar[:]=MODEL
                ncvar = ncOUT.createVariable('float','f',('depth', ))
                ncvar[:]=REF
                ncvar = ncOUT.createVariable('flag' ,'i',('depth', ))
                ncvar[:]=FLAG

                setattr(ncOUT, 'longitude', p.lon)
                setattr(ncOUT, 'latitude' , p.lat)
                ncOUT.close()
        CR = checkreport(line, nP,nexcl, flag, np.nan)
        return CR

    def phytoC_check(self, model, ref, depth, p):
        bad=(np.isnan(model)) | (np.isnan(ref))
        MODEL = model[~bad]
        REF   = ref [~bad]
        DEPTH= depth[~bad]
        nP = len(DEPTH)

        mydiff = np.abs(MODEL-REF)
        flag1_array = (mydiff > 5)  & (DEPTH<=200)
        flag1 = flag1_array.sum() > 5
        nexcl = 0

        line=""
        flag=np.nan
        if ( flag1) :#  |  flag2): 
            if flag1: flag=1
            # if flag2: flag=2
            # if (flag1  & flag2) : flag=3
            line="%s\t%s\t%s\t%s\t%s\n" %( p._my_float.wmo, p.time.strftime("%Y%m%d"), p.lon, p.lat , flag)
            nexcl=nP
            FLAG = np.zeros((nP), np.int)
            FLAG[flag1_array] =1
            if self.verboselevel ==1 :
                outncfile="%s%s_%s.nc"  %(self.outdir + 'P_c.', p.time.strftime("%Y%m%d"), p._my_float.wmo )
                ncOUT = netCDF4.Dataset(outncfile,'w')
                ncOUT.createDimension('depth',nP)

                ncvar = ncOUT.createVariable('model','f',('depth', ))
                ncvar[:]=MODEL
                ncvar = ncOUT.createVariable('float','f',('depth', ))
                ncvar[:]=REF
                ncvar = ncOUT.createVariable('flag' ,'i',('depth', ))
                ncvar[:]=FLAG

                setattr(ncOUT, 'longitude', p.lon)
                setattr(ncOUT, 'latitude' , p.lat)
                ncOUT.close()
        CR = checkreport(line, nP,nexcl, flag, np.nan)
        return CR


    def oxygen_check(self, model, ref, depth, p):
        bad=(np.isnan(model)) | (np.isnan(ref))
        MODEL = model[~bad]
        REF   = ref [~bad]
        DEPTH= depth[~bad]
        nP = len(DEPTH)

        REF_sat = oxy_sat(p)
 
        mydiff = np.abs(REF[0]-REF_sat)
        flag1 = (mydiff > 20)
        nexcl = 0

        line=""
        flag=np.nan
        if ( flag1) :
            if flag1: flag=1
            line="%s\t%s\t%s\t%s\t%s\n" %( p._my_float.wmo, p.time.strftime("%Y%m%d"), p.lon, p.lat , flag)
            nexcl=nP
            FLAG = np.zeros((1), np.int)
            FLAG[0] =1
            if self.verboselevel ==1 :
                outncfile="%s%s_%s.nc"  %(self.outdir + 'O2o.', p.time.strftime("%Y%m%d"), p._my_float.wmo )
                ncOUT = netCDF4.Dataset(outncfile,'w')
                ncOUT.createDimension('depth',nP)

                ncvar = ncOUT.createVariable('model','f',('depth', ))
                ncvar[:]=MODEL
                ncvar = ncOUT.createVariable('float','f',('depth', ))
                ncvar[:]=REF
                ncvar = ncOUT.createVariable('flag' ,'i',('depth', ))
                ncvar[:]=FLAG

                setattr(ncOUT, 'longitude', p.lon)
                setattr(ncOUT, 'latitude' , p.lat)
                ncOUT.close()
        CR = checkreport(line, nP,nexcl, flag, np.nan)
        return CR


    def perform(self, varname, Pres, Value, Qc, Model, p):
        ''' Works of Ref profiles '''
        
        CR=None
        if varname=='N3n' :
            CR = self.nitrate_check(Model, Value, Pres,  p)
            if CR.line != '':
#                 BAD__COUNTER[wmo] += 1
#                 OBS_BAD__COUNTER[wmo] += nexc
#                 OBS_GOOD_COUNTER[wmo] += nobs-nexc
                if CR.reason in [1,3,-1]:
                    Pres  = np.array([],np.float32)
                    Value = np.array([],np.float32)
                    Qc    = np.array([],np.int)
                if (CR.reason==2):
                    for iid in range(len(Pres)):
                        if Pres[iid]>CR.depthexc:
                            indexexc = iid
                            break
                    Pres =  Pres[:indexexc+1]
                    Value = Value[:indexexc+1]
                    Qc    = Qc[:indexexc+1]
#            else:
#                 GOOD_COUNTER[wmo] += 1
#                 OBS_GOOD_COUNTER[wmo] += Prof.shape[0]
        if varname=='P_l' :
            CR = self.chlorophyll_check(Model, Value, Pres, p)
            if CR.line != '':
                Pres  = np.array([],np.float32)
                Value = np.array([],np.float32)
                Qc    = np.array([],np.int)

        if varname=='P_c':
            CR = self.phytoC_check(Model, Value, Pres, p)
            if CR.line != '':
                Pres  = np.array([],np.float32)
                Value = np.array([],np.float32)
                Qc    = np.array([],np.int)

        if varname=='O2o':
            CR = self.oxygen_check(Model, Value, Pres, p)
            if CR.line != '':
                Pres  = np.array([],np.float32)
                Value = np.array([],np.float32)
                Qc    = np.array([],np.int)

        return Pres, Value, Qc, CR


if __name__=="__main__":
    import superfloat
    OUTDIR="/gpfs/scratch/userexternal/gbolzon0/ateruzzi/NITRATE_CHECK/nitrate_profiles/"
    A = check(OUTDIR,verboselevel=1)
    BF=superfloat.BioFloat.from_file('/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V5C/SUPERFLOAT/6901764/MR6901764_114.nc')
    p = superfloat.BioFloatProfile(BF.time,BF.lon, BF.lat, BF, BF.available_params,None)
    Pres, Value,Qc = p.read('NITRATE')
    Model = Value-10.3
    Pres,Values,Qc,CR = A.perform("N3n", Pres, Value, Qc, Model, p)
    #Pres,Values,Qc,CR = A.perform("P_l", Pres, Value, Qc, Model, p)

    

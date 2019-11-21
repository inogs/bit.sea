import numpy as np
import netCDF4

class check():
    def __init__(self,OUTDIR):
        self.outdir=OUTDIR

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
        
            outncfile="%s%s_%s.nc"  %(self.outdir + '/OUTNC/N3n.', p.time.strftime("%Y%m%d"), p._my_float.wmo )
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
        return line,firstdepth

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
           
    
        line=""
        if ( flag1) :#  |  flag2): 
            if flag1: flag=1
            # if flag2: flag=2
            # if (flag1  & flag2) : flag=3
            line="%s\t%s\t%s\t%s\t%s\n" %( p._my_float.wmo, p.time.strftime("%Y%m%d"), p.lon, p.lat , flag)
        
            FLAG = np.zeros((nP), np.int)
            FLAG[flag1_array] =1
            # FLAG[flag2_array] =2
        
            outncfile="%s%s_%s.nc"  %(OUTDIR + '/OUTNC/P_l.', p.time.strftime("%Y%m%d"), p._my_float.wmo )
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
        return line
    def perform(self, varname, Pres, Value, Qc, Model, p):
        ''' Works of Ref profiles '''
        
        if varname=='N3n' :
            line,depthexc = check.nitrate_check(Model, Value, Pres,  p)
            if line != '':
                #LINES.append(line)
                reason = np.int(line.rsplit('\t')[-1].rsplit('\n')[0])
#                nobs = np.int(line.rsplit('\t')[-3])
#                nexc = np.int(line.rsplit('\t')[-2])
#                 BAD__COUNTER[wmo] += 1
#                 OBS_BAD__COUNTER[wmo] += nexc
#                 OBS_GOOD_COUNTER[wmo] += nobs-nexc
                if reason in [1,3-1]:
                    Pres  = np.array([],np.float32)
                    Value = np.array([],np.float32)
                    Qc    = np.array([],np.int)
                if (reason==2):
                    for iid in range(len(Pres)):
                        if Pres[iid]>depthexc:
                            indexexc = iid
                            break
                    Pres =  Pres[:indexexc+1]
                    Value = Value[:indexexc+1]
                    Qc    = Qc[:indexexc+1]
#            else:
#                 GOOD_COUNTER[wmo] += 1
#                 OBS_GOOD_COUNTER[wmo] += Prof.shape[0]        
        return Pres, Value, Qc, line


if __name__=="__main__":
    OUTDIR="/gpfs/scratch/userexternal/gbolzon0/ateruzzi/NITRATE_CHECK/nitrate_profiles/"
    A = check(OUTDIR)
    
    
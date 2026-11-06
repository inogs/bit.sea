import numpy as np
import bitsea.Sat.SatManager as Sat
import bitsea.matchup.matchup as matchup
from bitsea.commons.dataextractor import DataExtractor
from bitsea.basins import V2 as OGS #MODIFICARE CON I SOTTOBACINI del NAd (al momento è adr1)
import scipy.io.netcdf as NC
import os

def weighted_mean(Conc, Weight):

    Weight_sum      = Weight.sum()
    Mass            = (Conc * Weight).sum()
    Weighted_Mean   = Mass/Weight_sum
    return Weighted_Mean

def weighted_var(Conc, Weight):
    ConcMean = weighted_mean(Conc,Weight)
    Weight_sum      = Weight.sum()
    Mass            = ((Conc-ConcMean)**2 * Weight).sum()
    Weighted_Var    = Mass/Weight_sum
    return Weighted_Var

def SatValidation(var, modfile,satfile,climafile,TheMask,outfile,SUB,COASTNESS_LIST,COASTNESS,nSUB):
    assert var in ['P_l','kd490']
    varname_sat={'P_l':'CHL','kd490':'KD490'}
    if not os.path.exists(satfile):
        print("Not existing file: " + satfile)
        return
    
    print(outfile)
    nCOAST = len(COASTNESS_LIST)
    De = DataExtractor(TheMask,filename=modfile,varname=var,dimvar=2)
    Model = De.values[:,:]

    if climafile is not None:
        Declim = DataExtractor(TheMask,filename=climafile,varname=var,dimvar=2)
        Clima = Declim.values[:,:]
#    try:
#        Sat24 = Sat.readfromfile(satfile,varname_sat[var]) # weekly
#    except:
        Sat24 = Sat.convertinV4format(Sat.readfromfile(satfile, varname_sat[var]))  # daily

    cloudsLand = np.isnan(Sat24)
    Sat24[cloudsLand] = -999.0
    cloudsLand = Sat24==-999.0
    modelLand  = Model > 1.0e+19
    nodata     = cloudsLand | modelLand



    VALID_POINTS                       = np.zeros((nSUB,nCOAST),np.float32)

    MODEL_MEAN                         = np.zeros((nSUB,nCOAST),np.float32)
    SAT___MEAN                         = np.zeros((nSUB,nCOAST),np.float32)
    MODEL_VARIANCE                     = np.zeros((nSUB,nCOAST),np.float32)
    SAT___VARIANCE                     = np.zeros((nSUB,nCOAST),np.float32)
    BGC_CLASS4_CHL_BIAS_SURF_BASIN     = np.zeros((nSUB,nCOAST),np.float32)
    BGC_CLASS4_CHL_RMS_SURF_BASIN      = np.zeros((nSUB,nCOAST),np.float32)
    BGC_CLASS4_CHL_CORR_SURF_BASIN     = np.zeros((nSUB,nCOAST),np.float32)


    ANOMALY_CORR     = np.zeros((nSUB,nCOAST),np.float32)

    for icoast, coast in enumerate(COASTNESS_LIST):

#        for isub, sub in enumerate(OGS.P):
        for isub, sub in enumerate(OGS.adr):
            selection = SUB[sub.name] & (~nodata) & COASTNESS[coast]
            M = matchup.matchup(Model[selection], Sat24[selection])
            VALID_POINTS[isub,icoast] = M.number()
            if M.number() > 0 :
                BGC_CLASS4_CHL_RMS_SURF_BASIN[isub,icoast]  = M.RMSE()
                BGC_CLASS4_CHL_BIAS_SURF_BASIN[isub,icoast] = M.bias()
                BGC_CLASS4_CHL_CORR_SURF_BASIN[isub,icoast] = M.correlation()

                weight = TheMask.area[selection]
                MODEL_MEAN[isub,icoast] = weighted_mean( M.Model,weight)
                SAT___MEAN[isub,icoast] = weighted_mean( M.Ref,  weight)
                MODEL_VARIANCE[isub,icoast] = weighted_var( M.Model,weight)
                SAT___VARIANCE[isub,icoast] = weighted_var( M.Ref,  weight)



                if climafile is not None:
                    Mclima = matchup.matchup(Model[selection]-Clima[selection],  \
                                            Sat24[selection]-Clima[selection])
                    ANOMALY_CORR[isub,icoast] = Mclima.correlation()

            else:
                BGC_CLASS4_CHL_RMS_SURF_BASIN[isub,icoast] = np.nan
                BGC_CLASS4_CHL_BIAS_SURF_BASIN[isub,icoast]= np.nan
                BGC_CLASS4_CHL_CORR_SURF_BASIN[isub,icoast]= np.nan
                MODEL_MEAN[isub,icoast] = np.nan
                SAT___MEAN[isub,icoast] = np.nan
                MODEL_VARIANCE[isub,icoast] = np.nan
                SAT___VARIANCE[isub,icoast] = np.nan

                ANOMALY_CORR[isub,icoast]     = np.nan


    ncOUT = NC.netcdf_file(outfile,'w')
    ncOUT.createDimension('nsub', nSUB)
    ncOUT.createDimension('ncoast', nCOAST)
    s=''
    for sub in OGS.P: s =s+sub.name + ","
    setattr(ncOUT,'sublist',s)
    s=''
    for coast in COASTNESS_LIST: s =s+coast + ","
    setattr(ncOUT,'coastlist',s)


    ncvar= ncOUT.createVariable('number', 'i', ('nsub','ncoast'))
    ncvar[:]= VALID_POINTS

    ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_RMS_SURF_BASIN', 'f', ('nsub','ncoast'))
    ncvar[:] = BGC_CLASS4_CHL_RMS_SURF_BASIN

    ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_BIAS_SURF_BASIN', 'f', ('nsub','ncoast'))
    ncvar[:] = BGC_CLASS4_CHL_BIAS_SURF_BASIN

    ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_CORR_SURF_BASIN', 'f', ('nsub','ncoast'))
    ncvar[:] = BGC_CLASS4_CHL_CORR_SURF_BASIN

    ncvar = ncOUT.createVariable('MODEL_MEAN', 'f', ('nsub','ncoast'))
    ncvar[:] = MODEL_MEAN

    ncvar=ncOUT.createVariable('SAT___MEAN', 'f', ('nsub','ncoast'))
    ncvar[:] = SAT___MEAN

    ncvar = ncOUT.createVariable('MODEL_VARIANCE', 'f', ('nsub','ncoast'))
    ncvar[:] = MODEL_VARIANCE

    ncvar=ncOUT.createVariable('SAT___VARIANCE', 'f', ('nsub','ncoast'))
    ncvar[:] = SAT___VARIANCE


    if climafile is not None:
        ncvar=ncOUT.createVariable('ANOMALY_CORRELATION', 'f', ('nsub','ncoast'))
        ncvar[:] = ANOMALY_CORR

    ncOUT.close()


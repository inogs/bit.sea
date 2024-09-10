from bitsea.instruments import bio_float
from bitsea.Float import oxygen_saturation
from bitsea.commons.time_interval import TimeInterval
from bitsea.basins.region import Rectangle
from bitsea.commons.utils import addsep
import os
import scipy.io.netcdf as NC
import numpy as np
import seawater as sw
import datetime

def remove_bad_sensors(Profilelist,var):
    '''

    Subsetter, filtering out bad sensors for that var

     Arguments:
      * Profilelist * list of Profile objects
      * var         * string

      Returns:
        a list of Profile Objects
    '''

    OUT_N3n = ["6903197","6901767","6901773","6901771"]
    OUT_O2o = ["6901510"]
    OUT_O2o = ["6901766",'6903235','6902902',"6902700"]
    # 0 6901766 has negative values

    if ( var == 'SR_NO3' ):
        return [p for p in Profilelist if p.name() not in OUT_N3n]

    if ( var == 'DOXY' ):
        return [p for p in Profilelist if p.name() not in OUT_O2o]

    return Profilelist

def dump_oxygen_file(outfile, p, Pres, Value, Qc, metadata, mode='a'):
    nP=len(Pres)
    if mode=='a':
        command = "cp %s %s.tmp" %(outfile,outfile)
        os.system(command)
    ncOUT = NC.netcdf_file(outfile + ".tmp" ,mode)

    if mode=='w': # if not existing file, we'll put header, TEMP, PSAL
        setattr(ncOUT, 'origin'     , 'coriolis')
        setattr(ncOUT, 'file_origin', metadata.filename)
        PresT, Temp, QcT = p.read('TEMP')  #, read_adjusted=False)
        PresT, Sali, QcS = p.read('PSAL')  # , read_adjusted=False)
        ncOUT.createDimension("DATETIME",14)
        ncOUT.createDimension("NPROF", 1)
        ncOUT.createDimension('nTEMP', len(PresT))
        ncOUT.createDimension('nPSAL', len(PresT))

        ncvar=ncOUT.createVariable("REFERENCE_DATE_TIME", 'c', ("DATETIME",))
        ncvar[:]=p.time.strftime("%Y%m%d%H%M%S")
        ncvar=ncOUT.createVariable("JULD", 'd', ("NPROF",))
        ncvar[:]=0.0
        ncvar=ncOUT.createVariable("LONGITUDE", "d", ("NPROF",))
        ncvar[:] = p.lon.astype(np.float64)
        ncvar=ncOUT.createVariable("LATITUDE", "d", ("NPROF",))
        ncvar[:] = p.lat.astype(np.float64)

        ncvar=ncOUT.createVariable('TEMP','f',('nTEMP',))
        ncvar[:]=Temp
        setattr(ncvar, 'variable'   , 'TEMP')
        setattr(ncvar, 'units'      , "degree_Celsius")
        ncvar=ncOUT.createVariable('PRES_TEMP','f',('nTEMP',))
        ncvar[:]=PresT
        ncvar=ncOUT.createVariable('TEMP_QC','f',('nTEMP',))
        ncvar[:]=QcT

        ncvar=ncOUT.createVariable('PSAL','f',('nTEMP',))
        ncvar[:]=Sali
        setattr(ncvar, 'variable'   , 'SALI')
        setattr(ncvar, 'units'      , "PSS78")
        ncvar=ncOUT.createVariable('PRES_PSAL','f',('nTEMP',))
        ncvar[:]=PresT
        ncvar=ncOUT.createVariable('PSAL_QC','f',('nTEMP',))
        ncvar[:]=QcS


    #print ("dumping oxygen on " + outfile)
    doxy_already_existing="nDOXY" in ncOUT.dimensions.keys()
    if not doxy_already_existing : ncOUT.createDimension('nDOXY', nP)
    ncvar=ncOUT.createVariable("PRES_DOXY", 'f', ('nDOXY',))
    ncvar[:]=Pres
    ncvar=ncOUT.createVariable("DOXY", 'f', ('nDOXY',))
    ncvar[:]=Value
    #if not doxy_already_existing:
    setattr(ncvar, 'status_var' , metadata.origin('DOXY').status_var)
    setattr(ncvar, 'variable'   , 'DOXY')
    setattr(ncvar, 'units'      , "mmol/m3")
    ncvar=ncOUT.createVariable("DOXY_QC", 'f', ('nDOXY',))
    ncvar[:]=Qc
    ncOUT.close()

    os.system("mv " + outfile + ".tmp " + outfile)


def get_outfile(p,outdir):
    wmo=p._my_float.wmo
    filename="%s%s/%s" %(outdir,wmo, os.path.basename(p._my_float.filename))
    return filename


def doxy_algorithm(p, outfile, metadata,writing_mode):
    Pres, Value, Qc = read_doxy(p)
    if Pres is None: return


    if p._my_float.status_var('DOXY')=='D':
        metadata.status_var='D'
        condition_to_write = True
    else:
        metadata.status_var='A'
        condition_to_write =  oxygen_saturation.oxy_check(Pres,Value,p)


    if condition_to_write:
        os.system('mkdir -p ' + os.path.dirname(outfile))
        dump_oxygen_file(outfile, p, Pres, Value, Qc, metadata,mode=writing_mode)
    else:
        print "Saturation Test not passed"

def read_doxy(pCor):
    Pres, Value, Qc = pCor.read('DOXY',read_adjusted=True)
    nP=len(Pres)
    if nP<5 :
        print "few values for " + pCor._my_float.filename
        return None, None, None
    ValueCconv=convert_oxygen(pCor, Pres, Value)
    return Pres, ValueCconv, Qc


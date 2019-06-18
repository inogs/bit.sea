import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates superfloat files of chla.
    Reads from Coriolis and LOV datasets.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--datestart','-s',
                                type = str,
                                required = True,
                                help = '''date in "%Y%m%d" format, e.g .20120101  ''')
    parser.add_argument(   '--dateend','-e',
                                type = str,
                                required = True,
                                help = '''date in "%Y%m%d" format , e.g 20200101 ''')
    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                default = "/gpfs/scratch/userexternal/gbolzon0/SuperFloat/",
                                help = 'path of the Superfloat dataset ')

    return parser.parse_args()

args = argument()



from instruments import bio_float
from instruments import lovbio_float
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import superfloat_generator
from commons.utils import addsep
import os,sys
import scipy.io.netcdf as NC
import numpy as np



def get_info(p,outdir):
    wmo=p._my_float.wmo
    filename="%s%s/MR%s_%03d.nc" %(outdir,wmo, wmo,p._my_float.cycle)
    return filename

def dumpfile(outfile,p_pos, p,Pres,chl_profile,Qc,metadata):
    PresT, Temp, QcT = p.read('TEMP', read_adjusted=False)
    PresT, Sali, QcS = p.read('PSAL', read_adjusted=False)

    print "dumping chla on " + outfile
    ncOUT = NC.netcdf_file(outfile,"w")
    ncOUT.createDimension("DATETIME",14)
    ncOUT.createDimension("NPROF", 1)
    ncOUT.createDimension('nTEMP', len(PresT))
    ncOUT.createDimension('nPSAL', len(PresT))
    ncOUT.createDimension('nCHLA', len(Pres ))    
    
    ncvar=ncOUT.createVariable("REFERENCE_DATE_TIME", 'c', ("DATETIME",))
    ncvar[:]=p_pos.time.strftime("%Y%m%d%H%M%S")
    ncvar=ncOUT.createVariable("JULD", 'd', ("NPROF",))
    ncvar[:]=0.0
    ncvar=ncOUT.createVariable("LONGITUDE", "d", ("NPROF",))
    ncvar[:] = p_pos.lon.astype(np.float64)
    ncvar=ncOUT.createVariable("LATITUDE", "d", ("NPROF",))
    ncvar[:] = p_pos.lat.astype(np.float64)
 
    
    ncvar=ncOUT.createVariable('TEMP','f',('nTEMP',))
    ncvar[:]=Temp
    setattr(ncvar, 'origin'     , metadata.origin)
    setattr(ncvar, 'file_origin', metadata.filename)
    setattr(ncvar, 'variable'   , 'TEMP')
    setattr(ncvar, 'units'      , "degree_Celsius")

    ncvar=ncOUT.createVariable('PRES_TEMP','f',('nTEMP',))
    ncvar[:]=PresT
    ncvar=ncOUT.createVariable('TEMP_QC','f',('nTEMP',))
    ncvar[:]=QcT
    
    ncvar=ncOUT.createVariable('PSAL','f',('nTEMP',))
    ncvar[:]=Sali
    setattr(ncvar, 'origin'     , metadata.origin)
    setattr(ncvar, 'file_origin', metadata.filename)
    setattr(ncvar, 'variable'   , 'SALI')
    setattr(ncvar, 'units'      , "PSS78")

    ncvar=ncOUT.createVariable('PRES_PSAL','f',('nTEMP',))
    ncvar[:]=PresT
    ncvar=ncOUT.createVariable('PSAL_QC','f',('nTEMP',))
    ncvar[:]=QcS

    ncvar=ncOUT.createVariable('CHLA','f',('nCHLA',))
    ncvar[:]=chl_profile
    setattr(ncvar, 'origin'     , metadata.origin)
    setattr(ncvar, 'file_origin', metadata.filename)
    setattr(ncvar, 'variable'   , 'CHLA_ADJUSTED')
    setattr(ncvar, 'units'      , "milligram/m3")

    ncvar=ncOUT.createVariable('PRES_CHLA','f',('nCHLA',))
    ncvar[:]=Pres
    ncvar=ncOUT.createVariable('CHLA_QC','f',('nCHLA',))
    ncvar[:]=Qc
    ncOUT.close()

def is_only_lov(pCor):
    if pCor is None: return True
    return not superfloat_generator.exist_variable('CHLA', pCor._my_float.filename)


OUTDIR = addsep(args.outdir)
TI     = TimeInterval(args.datestart,args.dateend,'%Y%m%d')
R = Rectangle(-6,36,30,46)

PROFILES_LOV =lovbio_float.FloatSelector('CHLA', TI, R)
wmo_list= lovbio_float.get_wmo_list(PROFILES_LOV)



for wmo in wmo_list:
    print wmo
    Profilelist = lovbio_float.filter_by_wmo(PROFILES_LOV, wmo)
    for ip, pLov in enumerate(Profilelist):
        pCor = bio_float.from_lov_profile(pLov, verbose=True)
        is_only_LOV = is_only_lov(pCor)
        if is_only_LOV:
            outfile = get_info(pLov,OUTDIR)
        else:
            outfile = get_info(pCor,OUTDIR)
        if os.path.exists(outfile): continue
        os.system('mkdir -p ' + os.path.dirname(outfile))

        if is_only_lov:
            Pres, Profile, Qc = pLov.read('CHLA',read_adjusted=True)
            profile_for_data = pLov
            if pCor is None:
                profile_for_pos= pLov
            else:
                profile_for_pos= pCor
            metadata = superfloat_generator.Metadata('lov', pLov._my_float.filename)

        else:
            #Pres, Profile,Qc = superfloat_generator.synthesis_profile(pLov, pCor)
            Pres, Profile,Qc = superfloat_generator.treating_coriolis(pCor)
            profile_for_data = pCor
            profile_for_pos  = pCor
            metadata = superfloat_generator.Metadata('Coriolis', pCor._my_float.filename)


        if Pres is None: continue # no data

        dumpfile(outfile, profile_for_pos, profile_for_data, Pres, Profile, Qc, metadata)


print "**************** Profiles only Coriolis *******************"
# Deve gestire i profili only_coriolis
PROFILES_COR =bio_float.FloatSelector('CHLA', TI, R)
wmo_list= bio_float.get_wmo_list(PROFILES_COR)

for wmo in wmo_list:
    print wmo
    Profilelist = bio_float.filter_by_wmo(PROFILES_COR, wmo)
    for ip, pCor in enumerate(Profilelist):
        outfile = get_info(pCor, OUTDIR)
        if os.path.exists(outfile): continue
        os.system('mkdir -p ' + os.path.dirname(outfile))
        Pres, CHL, Qc= superfloat_generator.treating_coriolis(pCor)
        metadata = superfloat_generator.Metadata('Coriolis', pCor._my_float.filename)
        if Pres is None: continue # no data
        dumpfile(outfile, pCor, pCor, Pres, CHL, Qc, metadata)


# DRIFTS in 6900807 7900591
    

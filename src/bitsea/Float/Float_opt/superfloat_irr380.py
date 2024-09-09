import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates superfloat files of irr_380.
    Reads from Float_opt_2019.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--datestart','-s',
                                type = str,
                                required = True,
                                help = '''date in yyyymmdd format''')
    parser.add_argument(   '--dateend','-e',
                                type = str,
                                required = True,
                                help = '''date in yyyymmdd format ''')
    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                default = "/gpfs/scratch/userexternal/gbolzon0/SUPERFLOAT/",
                                help = 'path of the Superfloat dataset ')
    parser.add_argument(   '--force', '-f',
                                action='store_true',
                                help = """Overwrite existing files
                                """)

    return parser.parse_args()

args = argument()

from instruments import optbio_float_2019
from instruments import optbio_float_2020
from commons.time_interval import TimeInterval
from basins.region import Rectangle
from Float import superfloat_generator
from commons.utils import addsep
import os
import scipy.io.netcdf as NC
import numpy as np
import seawater as sw


def dump_irr_380_file(outfile, p, Pres, Value, Qc, metatata, mode='w'):
    nP=len(Pres)
    if mode=='a':
        command = "cp %s %s.tmp" %(outfile,outfile)
        os.system(command)
    ncOUT = NC.netcdf_file(outfile + ".tmp" ,mode)

    if mode=='w': # if not existing file, we'll put header, TEMP, PSAL
        PresT, Temp, QcT = p.read('TEMP', read_adjusted=False)
        PresT, Sali, QcS = p.read('PSAL', read_adjusted=False)        
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

    print "dumping irr_380 on " + outfile
    irr_380_already_existing="nIRR_380" in ncOUT.dimensions.keys()
    if not irr_380_already_existing : ncOUT.createDimension('nIRR_380', nP)
    ncvar=ncOUT.createVariable("PRES_IRR_380", 'f', ('nIRR_380',))
    ncvar[:]=Pres
    ncvar=ncOUT.createVariable("IRR_380", 'f', ('nIRR_380',))
    ncvar[:]=Value
    if not irr_380_already_existing:
        setattr(ncvar, 'origin'     , metadata.origin)
        setattr(ncvar, 'file_origin', metadata.filename)
        setattr(ncvar, 'variable'   , 'IRR_380')
        setattr(ncvar, 'units'      , "W/m^2/nm")
    ncvar=ncOUT.createVariable("IRR_380_QC", 'f', ('nIRR_380',))
    ncvar[:]=Qc
    ncOUT.close()

    os.system("mv " + outfile + ".tmp " + outfile)


OUTDIR = addsep(args.outdir)
TI     = TimeInterval(args.datestart,args.dateend,'%Y%m%d')
R = Rectangle(-6,36,30,46)
force_writing_irr_380=args.force

PROFILES_OPT_2019 =optbio_float_2019.FloatSelector('IRR_380', TI, R)

wmo_list= optbio_float_2019.get_wmo_list(PROFILES_OPT_2019)
nWMOS=len(wmo_list)
#for p in PROFILES_OPT_2019:
#    p_2020=optbio_float_2020.from_profile(p, verbose=True)
#    import sys
#    sys.exit()


def get_outfile(p,outdir):
    filename="%s%s" %(outdir,p._my_float.filename)
    return filename

for iwmo, wmo in enumerate(wmo_list):
    print wmo, iwmo, " of ", nWMOS
    Profilelist=optbio_float_2019.filter_by_wmo(PROFILES_OPT_2019, wmo)
    for ip, p in enumerate(Profilelist):
        p_2020 = optbio_float_2020.from_profile(p, verbose=True)
        if p_2020 is None: continue
        outfile = get_outfile(p_2020,OUTDIR)
        metadata = superfloat_generator.Metadata('float_opt_2019', p._my_float.filename)
        #os.system('mkdir -p ' + os.path.dirname(outfile))

        if superfloat_generator.exist_valid(outfile):
            if not superfloat_generator.exist_variable('IRR_380', outfile):
                Pres, Value, Qc = p.read('IRR_380')
                if Pres is not None: dump_irr_380_file(outfile, p, Pres, Value, Qc, metadata,mode='a')
            else:
                if force_writing_irr_380:
                    Pres, Value, Qc = p.read('IRR_380')
                    if Pres is not None: dump_irr_380_file(outfile, p, Pres, Value, Qc, metadata,mode='a')
        else:
            print outfile + "not existing"
             #Pres, Value, Qc = p.read('RR_380')
             #if Pres is not None: dump_irr_380_file(outfile, p, Pres, Value, Qc, metadata,mode='w')


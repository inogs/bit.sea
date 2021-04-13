import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates superfloat files of nitrate.
    Reads from CORIOLIS dataset.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--datestart','-s',
                                type = str,
                                required = False,
                                help = '''date in yyyymmdd format''')
    parser.add_argument(   '--dateend','-e',
                                type = str,
                                required = False,
                                help = '''date in yyyymmdd format''')
    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                default = "/gpfs/scratch/userexternal/gbolzon0/SUPERFLOAT/",
                                help = 'path of the Superfloat dataset ')
    parser.add_argument(   '--force', '-f',
                                action='store_true',
                                help = """Overwrite existing files
                                """)
    parser.add_argument(   '--update_file','-u',
                                type = str,
                                required = False,
                                default = 'NO_file',
                                help = '''file with updated floats''')

    return parser.parse_args()

args = argument()

if (args.datestart == 'NO_data') & (args.dateend == 'NO_data') & (args.update_file == 'NO_file'):
    raise ValueError("No file nor data inserted: you have to pass either datastart and dataeend or the update_file")

if ((args.datestart == 'NO_data') or (args.dateend == 'NO_data')) & (args.update_file == 'NO_file'):
    raise ValueError("No file nor data inserted: you have to pass both datastart and dataeend")

from instruments import bio_float
from Float.canyon_b_N3n import canyon_nitrate_correction
from Float.woa_N3n import woa_nitrate_correction
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import superfloat_generator
from commons.utils import addsep
import os
import scipy.io.netcdf as NC
import numpy as np
from commons.layer import Layer
import seawater as sw
import datetime

import basins.V2 as basV2
from static.climatology import get_climatology
SUBLIST = basV2.P.basin_list
LayerList=[Layer(400,600),Layer(600,800),Layer(800,1000)]
N3n_clim, N3n_std = get_climatology('N3n', SUBLIST, LayerList)


def get_outfile(p,outdir):
    wmo=p._my_float.wmo
    filename="%s%s/%s" %(outdir,wmo, os.path.basename(p._my_float.filename))
    return filename

def convert_nitrate(p,pres,profile):
    '''
    from micromol/Kg to  mmol/m3
    '''
    if pres.size == 0: return profile
    Pres, temp, Qc = p.read("TEMP",read_adjusted=False)
    Pres, sali, Qc = p.read("PSAL",read_adjusted=False)
    density = sw.dens(sali,temp,Pres)
    density_on_z = np.interp(pres,Pres,density)
    return profile * density_on_z/1000.

def dump_nitrate_file(outfile, p, Pres, Value, Qc, metadata,mode='w'):
    
    nP=len(Pres)
    if mode=='a':
        command = "cp %s %s.tmp" %(outfile,outfile)
        os.system(command)
    ncOUT = NC.netcdf_file(outfile + ".tmp",mode)
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

    print "dumping nitrate on " + outfile
    nitrate_already_existing="nNITRATE" in ncOUT.dimensions.keys()
    if not nitrate_already_existing: ncOUT.createDimension('nNITRATE', nP)
    ncvar=ncOUT.createVariable("PRES_NITRATE", 'f', ('nNITRATE',))
    ncvar[:]=Pres
    ncvar=ncOUT.createVariable("NITRATE", 'f', ('nNITRATE',))
    ncvar[:]=Value
    if not nitrate_already_existing:
        setattr(ncvar, 'origin'     , metadata.origin)
        setattr(ncvar, 'file_origin', metadata.filename)
        setattr(ncvar, 'variable'   , 'SR_NO3_ADJUSTED')
        setattr(ncvar, 'units'      , "mmol/m3")
    ncvar=ncOUT.createVariable("NITRATE_QC", 'f', ('nNITRATE',))
    ncvar[:]=Qc
    ncOUT.close()
    os.system("mv " + outfile + ".tmp " + outfile)

OUTDIR = addsep(args.outdir)
force_writing_nitrate=args.force
input_file=args.update_file

if input_file == 'NO_file':
    TI = TimeInterval(args.datestart,args.dateend,'%Y%m%d')
    R = Rectangle(-6,36,30,46)


    PROFILES_COR =bio_float.FloatSelector('NITRATE', TI, R)
    wmo_list= bio_float.get_wmo_list(PROFILES_COR)


    for wmo in wmo_list:
        print wmo
        Profilelist = bio_float.filter_by_wmo(PROFILES_COR, wmo)
        for ip, pCor in enumerate(Profilelist):
            outfile = get_outfile(pCor,OUTDIR)
        F=pCor._my_float
        writing_mode='w'
        if superfloat_generator.exist_valid(outfile): writing_mode='a'

        condition_to_write = ~superfloat_generator.exist_valid_variable('NITRATE',outfile)
        if force_writing_nitrate: condition_to_write=True

        metatata = superfloat_generator.Metadata('coriolis',F.filename)

        if not condition_to_write: continue
        if pCor._my_float.status_var('NITRATE')=='R': continue

        Pres, Value, Qc= pCor.read("NITRATE", read_adjusted=True)
        nP=len(Pres)
        if nP<5 :
            print "few values for " + F.filename
            continue
        if Pres[-1]<100:
            print "depth < 100 for "+ F.filename
            continue

        os.system('mkdir -p ' + os.path.dirname(outfile))
        if pCor._my_float.status_var('NITRATE')=='D':
            dump_nitrate_file(outfile, pCor, Pres, Value, Qc, metatata,mode=writing_mode)


        if pCor._my_float.status_var('NITRATE')=='A':
            if superfloat_generator.exist_valid_variable('DOXY', outfile):
                DOXYp, DOXY, _ = pCor.read('DOXY',read_adjusted=True)


            Pres, Value, Qc, t_lev, nit = canyon_nitrate_correction(pCor, Pres, Value, Qc, DOXYp, DOXY)
            outOfClimatology = False
            for ilayer, layer in enumerate(LayerList):
                if (t_lev >= layer.top) & (t_lev < layer.bottom) :
                    for iSub,sub in enumerate(SUBLIST):
                        if sub.is_inside(pCor.lon,pCor.lat):
                            if np.abs(N3n_clim[iSub,ilayer] - nit ) > 2 : outOfClimatology = True

            if outOfClimatology:
                print "Out of climatology for " + F.filename
                continue
            dump_nitrate_file(outfile, pCor, Pres, Value, Qc, metatata,mode=writing_mode)
        else:
            if Pres[-1]>600:
                Pres, Values, Qc = woa_nitrate_correction(pCor)
                dump_nitrate_file(outfile, pCor, Pres, Value, Qc, metatata,mode=writing_mode)
            else:
                print "WOA correction not applicable for max(depth) < 600 m"

else:
    OUTDIR = addsep(args.outdir)
#    force_writing_oxygen=args.force
    mydtype= np.dtype([
        ('file_name','S200'),
        ('date','S200'),
        ('latitude',np.float32),
        ('longitude',np.float32),
        ('ocean','S10'),
        ('profiler_type',np.int),
        ('institution','S10'),
        ('parameters','S200'),
        ('parameter_data_mode','S100'),
        ('date_update','S200')] )

    INDEX_FILE=np.loadtxt(input_file,dtype=mydtype, delimiter=",",ndmin=1,skiprows=0)
    nFiles=INDEX_FILE.size

    for iFile in range(nFiles):
        timestr          = INDEX_FILE['date'][iFile]
        lon              = INDEX_FILE['longitude' ][iFile]
        lat              = INDEX_FILE['latitude' ][iFile]
        filename         = INDEX_FILE['file_name'][iFile]
        available_params = INDEX_FILE['parameters'][iFile]
        parameterdatamode= INDEX_FILE['parameter_data_mode'][iFile]
        float_time = datetime.datetime.strptime(timestr,'%Y%m%d%H%M%S')
        filename=filename.replace('coriolis/','').replace('profiles/','')


        if 'NITRATE' not in available_params: continue

        pCor=bio_float.profile_gen(lon, lat, float_time, filename, available_params,parameterdatamode)
        wmo=pCor._my_float.wmo

        outfile = get_outfile(pCor,OUTDIR)
        F=pCor._my_float
        writing_mode='w'
        if superfloat_generator.exist_valid(outfile): writing_mode='a'

        metatata = superfloat_generator.Metadata('coriolis',F.filename)

        if pCor._my_float.status_var('NITRATE')=='R': continue

        Pres, Value, Qc= pCor.read("NITRATE", read_adjusted=True)
        nP=len(Pres)
        if nP<5 :
            print "few values for " + F.filename
            continue
        if Pres[-1]<100:
            print "few values for " + F.filename
            continue

        os.system('mkdir -p ' + os.path.dirname(outfile))
        if pCor._my_float.status_var('NITRATE')=='D':
            dump_nitrate_file(outfile, pCor, Pres, Value, Qc, metatata,mode=writing_mode)


        if pCor._my_float.status_var('NITRATE')=='A':
            if superfloat_generator.exist_valid_variable('DOXY', outfile):
                DOXYp, DOXY, _ = pCor.read('DOXY',read_adjusted=True)


                Pres, Value, Qc, t_lev, nit = canyon_nitrate_correction(pCor, Pres, Value, Qc, DOXYp, DOXY)
                outOfClimatology = False
                for ilayer, layer in enumerate(LayerList):
                    if (t_lev >= layer.top) & (t_lev < layer.bottom) :
                        for iSub,sub in enumerate(SUBLIST):
                            if sub.is_inside(pCor.lon,pCor.lat):
                                if np.abs(N3n_clim[iSub,ilayer] - nit ) > 2 : outOfClimatology = True

                if outOfClimatology:
                        print "Out of climatology for " + F.filename
                        continue
                dump_nitrate_file(outfile, pCor, Pres, Value, Qc, metatata,mode=writing_mode)
            else:
                if Pres[-1]>600:
                    Pres, Values, Qc = woa_nitrate_correction(pCor)
                    dump_nitrate_file(outfile, pCor, Pres, Value, Qc, metatata,mode=writing_mode)
                else:
                    print "WOA correction not applicable for max(depth) < 600 m"



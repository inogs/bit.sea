import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates superfloat files of chla.
    Reads from CORIOLIS dataset.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--datestart','-s',
                                type = str,
                                required = True,
                                help = '''date in yyyymmdd format''')
    parser.add_argument(   '--dateend','-e',
                                type = str,
                                required = True,
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

    return parser.parse_args()

args = argument()

from instruments import bio_float
from canyon_b_N3n import canyon_nitrate_correction
from commons.time_interval import TimeInterval
from basins.region import Rectangle
import superfloat_generator
from commons.utils import addsep
import os
import scipy.io.netcdf as NC
import numpy as np
from commons.layer import Layer

import basins.V2 as basV2
from static.climatology import get_climatology
SUBLIST = basV2.P.basin_list
LayerList=[Layer(400,600),Layer(600,800),Layer(800,1000)]
N3n_clim, N3n_std = get_climatology('N3n', SUBLIST, LayerList)


def get_outfile(p,outdir):
    wmo=p._my_float.wmo
    filename="%s%s/MR%s_%03d.nc" %(outdir,wmo, wmo,p._my_float.cycle)
    return filename


def dump_nitrate_file(outfile, p_pos, p, Pres, Value, Qc, metadata,mode='w'):
    
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
TI     = TimeInterval(args.datestart,args.dateend,'%Y%m%d')
R = Rectangle(-6,36,30,46)


PROFILES_COR =bio_float.FloatSelector('NITRATE', TI, R)
wmo_list= bio_float.get_wmo_list(PROFILES_COR)

force_writing_nitrate=args.force

for wmo in wmo_list:
    print wmo
     # should be filtered 6901653, 6901655, 6901657 6901649, 6901764 6903197 6901605
    #1529 1513 1511 1863 1862 1860 1864 1776 1775 0807 0591
    Profilelist = bio_float.filter_by_wmo(PROFILES_COR, wmo)
    for ip, pCor in enumerate(Profilelist):
        outfile = get_outfile(pCor,OUTDIR)
        profile_for_position=pCor

        metatata = superfloat_generator.Metadata('coriolis',pCor._my_float.filename)

        F=pCor._my_float
        F.load_basics()
        if F.status_var('DOXY') == 'R':
            print "DOXY STATUS = R for", F.filename
            continue

        Pres, Value, Qc= pCor.read("NITRATE", read_adjusted=True)
        nP=len(Pres)
        if nP<5 :
            print "few values for " + F.filename
            continue
        if Pres[-1]<100:
            print "depth < 100 for "+ F.filename
            continue
        DOXYp, DOXY, _ = pCor.read('DOXY',True)
        nP=len(DOXYp)
        if nP < 5:
            print "Not enough DOXY_ADJUSTED values for " + F.filename
            continue
        os.system('mkdir -p ' + os.path.dirname(outfile))

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

        if superfloat_generator.exist_valid(outfile):
            if not superfloat_generator.exist_variable('NITRATE', outfile):
                dump_nitrate_file(outfile, profile_for_position, pCor, Pres, Value, Qc, metatata,mode='a')
            else:
                if force_writing_nitrate:
                    dump_nitrate_file(outfile, profile_for_position, pCor, Pres, Value, Qc, metatata,mode='a')
        else:
            print outfile + " not found"
            dump_nitrate_file(outfile, profile_for_position, pCor,Pres, Value, Qc, metatata,mode='w')


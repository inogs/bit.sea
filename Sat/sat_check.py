import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Apply check based on climatology to sat ORIG files.
    Produces CHECKED and REJECTED files for each date at satellite resolution.
    Statistics provided if argument dirstats is setted.
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--origdir', '-i',
                                type = str,
                                required = True,
                                help = ''' ORIG sat directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/ORIG/'''

                                )

    parser.add_argument(   '--checkdir', '-o',
                                type = str,
                                required = True,
                                help = ''' Base for CHECKED and REJECTED sat directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/'''

                                )

    parser.add_argument(   '--climfile', '-c',
                                type = str,
                                required = True,
                                help = ''' Climatology .nc file used to apply check on sat data, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/CCI02/SatClimatology.nc'''

                                )

    parser.add_argument(   '--mesh', '-m',
                                type = str,
                                required = True,
                                choices = ['SatOrigMesh','V4mesh','V1mesh','KD490mesh','SAT1km_mesh', 'Mesh24'],
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
                                )

    parser.add_argument(   '--submaskdir', '-s',
                                type = str,
                                required = True,
                                help = '''  Directory with subbasins on satellite mask'''
                                )

    parser.add_argument(   '--statsdir', '-w',
                                type = str,
                                required = False,
                                help = ''' STATS  directory'''

                                )

    return parser.parse_args()


args = argument()
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins import V2
from postproc import masks
import numpy as np
from commons.utils import addsep
import os
import scipy.io.netcdf as NC

import SatManager as Sat


ORIGDIR   = addsep(args.origdir)
SUBMASKDIR  = addsep(args.submaskdir)
CHECKDIR  = addsep(args.checkdir)
CLIM_FILE = args.climfile

if args.statsdir is not None:
    STATSDIR  = addsep(args.statsdir)

maskSat = getattr(masks,args.mesh)

reset = False

Timestart="19501231"
Time__end="20500101"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TL_orig = TimeList.fromfilenames(TI, ORIGDIR ,"*.nc",prefix='',dateformat='%Y%m%d')

somecheck = False
for iTime,filename in enumerate(TL_orig.filelist):
    outfile = CHECKDIR + '/CHECKED/' + os.path.basename(filename)
    rejfile = CHECKDIR + '/REJECTED/' + os.path.basename(filename)
    iDate = TL_orig.Timelist[iTime] 
    date8 = iDate.strftime('%Y%m%d')
    if (not os.path.exists(outfile)) or (not os.path.exists(rejfile)):
        somecheck = True
        break
    if args.statsdir is not None:
        for masktype in ['ORIG','CHECK']:
            filestats = STATSDIR + os.path.basename(filename)[:-3] + '_stats' + masktype + '.nc'
            if (not os.path.exists(filestats)):
                somecheck = True
                break
    if somecheck:
        break

if somecheck or reset:
    print('Read climatology')
    MEAN,STD = Sat.readClimatology(CLIM_FILE)

    filemaskmed = SUBMASKDIR + 'maskmed_1kmAll.npy'
    maskmed_1km = np.load(filemaskmed)

    masksub_M = {}
    nsub = 0
    subnames = ''
    for sub in V2.P:
        if sub.name=='atl': continue
        nsub += 1
        subnames += sub.name + ', '
        filemasksub = SUBMASKDIR + 'masksub.' + sub.name + 'All.npy'
        masksub = np.load(filemasksub)
        masksub_M[sub.name] = np.zeros((maskSat.jpj,maskSat.jpi),dtype=bool)
        masksub_M[sub.name][maskmed_1km] = masksub


else:
    print('All checks done')



for iTime, filename in enumerate(TL_orig.filelist):
    outfile = CHECKDIR + '/CHECKED/' + os.path.basename(filename)
    rejfile = CHECKDIR + '/REJECTED/' + os.path.basename(filename)
    exit_condition = os.path.exists(outfile) and os.path.exists(rejfile) and (not reset)
    exit_condmask = np.zeros(2,dtype=bool)
    iDate = TL_orig.Timelist[iTime] 
    date8 = iDate.strftime('%Y%m%d')
    if exit_condition:
        if args.statsdir is None:
            continue
        else:
            for imsk,masktype in enumerate(['ORIG','CHECK']):
                filestats = STATSDIR + os.path.basename(filename)[:-3] + '_stats' + masktype + '.nc'
                exit_condmask[imsk] = (os.path.exists(filestats)) and (not reset)
            if all(exit_condmask):
                continue

    julian = int( iDate.strftime("%j") )
    print(' ... day ' + np.str(julian) + '  of ' + np.str(iDate.year))
    if julian == 366:
        julian = 365


    maskreject = np.zeros_like(maskmed_1km,dtype=int)

    DAILY_REF_MEAN = MEAN[julian-1,:,:]
    DAILY_REF_STD  =  STD[julian-1,:,:]    
    
    CHL_IN = Sat.readfromfile(filename)

    if args.mesh == 'SatOrigMesh': CHL_IN[581:,164:] = Sat.fillValue # BLACK SEA

    cloudsLandTIME = CHL_IN         == Sat.fillValue
    cloudlandsCLIM = DAILY_REF_MEAN == Sat.fillValue
    
    CHL_OUT = CHL_IN.copy()
    CHL_OUT[cloudsLandTIME] = Sat.fillValue
    CHL_OUT[cloudlandsCLIM] = Sat.fillValue
    counter_refNAN = (~cloudsLandTIME & cloudlandsCLIM).sum(axis=None)
    maskreject[cloudsLandTIME] = -5
    maskreject[~cloudsLandTIME & cloudlandsCLIM] = 2
    
    outOfRange = np.abs(CHL_IN - DAILY_REF_MEAN) > DAILY_REF_STD *2.0
    outOfRange[cloudsLandTIME | cloudlandsCLIM ] = False
    
    counter_elim = outOfRange.sum(axis = None)
    CHL_OUT[outOfRange] = Sat.fillValue 
    maskreject[outOfRange] = 1
    
    if exit_condition==False:
        print 'Done check with ', filename, '  (',iTime+1,' of ', len(TL_orig.filelist), ')'
        print 'Rejection:  after check', counter_elim, ' values'
        print 'rejected for NAN in Climatology', counter_refNAN, ' values'
        Sat.dumpGenericNativefile(outfile, CHL_OUT, "CHL",mesh=maskSat)
        Sat.dumpGenericNativefile(rejfile, maskreject, "RejInd",mesh=maskSat)


    if (args.statsdir is not None) and (not all(exit_condmask)):
        masksubday = {}
        for sub in V2.P:
            if sub.name=='atl': continue
            masksubday[sub.name] = {}
            masksubday[sub.name]['ORIG'] = masksub_M[sub.name] & (cloudsLandTIME == False)
            masksubday[sub.name]['CHECK'] = masksub_M[sub.name] & (cloudsLandTIME == False) & (maskreject == 0)

        for masktype in ['ORIG','CHECK']:
            filestats = STATSDIR + os.path.basename(filename)[:-3] + '_stats' + masktype + '.nc'
            stats_clima = np.zeros((nsub,12))
            stats_clima[:,:] = np.nan
            stats_day = np.zeros((nsub,8))
            stats_day[:,:] = np.nan
            print(masktype + ' --- Cycle on sub   ---')
            for isub,sub in enumerate(V2.P):
                if sub.name=='atl': continue
                climadmean = DAILY_REF_MEAN[masksubday[sub.name][masktype]]
                climadmean[climadmean<0] = np.nan
                climadstd = DAILY_REF_STD[masksubday[sub.name][masktype]]
                climadstd[np.isnan(climadmean)] = np.nan

                climadmadd_2std = climadmean + 2*climadstd
                climadmsub_2std = climadmean - 2*climadstd

                chlsub = CHL_IN[masksubday[sub.name][masktype]]
                chlsub[chlsub<0] = np.nan

                if np.nansum(masksubday[sub.name][masktype])>0:
                    #print(sub.name + ' ... Np ' + np.str(np.nansum(masksubday)))
                    stats_clima[isub,0] = np.nanmean(climadmean)
                    stats_clima[isub,1] = np.nanmin(climadmean)
                    stats_clima[isub,2] = np.nanmax(climadmean)
                    stats_clima[isub,3] = np.nanmean(climadmadd_2std)
                    stats_clima[isub,4] = np.nanmin(climadmadd_2std)
                    stats_clima[isub,5] = np.nanmax(climadmadd_2std)
                    stats_clima[isub,6] = np.nanmean(climadmsub_2std)
                    stats_clima[isub,7] = np.nanmin(climadmsub_2std)
                    stats_clima[isub,8] = np.nanmax(climadmsub_2std)
                    stats_clima[isub,9] = np.nanmean(climadstd)
                    stats_clima[isub,10] = np.nanmin(climadstd)
                    stats_clima[isub,11] = np.nanmax(climadstd)


                    stats_day[isub,0] = np.nanmean(chlsub)
                    stats_day[isub,1] = np.nanmin(chlsub)
                    stats_day[isub,2] = np.nanmax(chlsub)
                    stats_day[isub,3] = np.nanmean(np.abs(chlsub-climadmean)/climadstd)
                    stats_day[isub,4] = np.nanmean(chlsub/climadmean)

                    stats_day[isub,5] = np.nansum(masksubday[sub.name]['ORIG'])
                    stats_day[isub,6] = np.nansum(outOfRange[masksubday[sub.name]['ORIG']])   #maskreject=1
                    stats_day[isub,7] = np.nansum(~cloudsLandTIME & cloudlandsCLIM & \
                                                        masksubday[sub.name]['ORIG']) #maskreject=2


            ncOUT = NC.netcdf_file(filestats,'w')
            ncOUT.createDimension('subbasin',nsub)
            ncOUT.createDimension('stattype_clim',12)
            ncOUT.createDimension('stattype_day',8)

            setattr(ncOUT,'sublist',subnames)

            ncvar = ncOUT.createVariable('SUBstatistics_clim','f',('subbasin','stattype_clim'))
            ncvar[:] = stats_clima

            ncvar = ncOUT.createVariable('SUBstatistics_day','f',('subbasin','stattype_day'))
            ncvar[:] = stats_day

            ncOUT.close()
            # fid = open(fileclim,'wb')
            # pickle.dump(stats_clima,fid)
            # fid.close()

            # fid = open(filechlsub,'wb')
            # pickle.dump(stats_day,fid)
            # fid.close()

        print 'Done statistics with ', filename, '  (',iTime+1,' of ', len(TL_orig.filelist), ')'
        print '   ---------------------------------------------------   '


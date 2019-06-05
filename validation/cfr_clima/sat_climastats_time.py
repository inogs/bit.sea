import numpy as np
import argparse
import pickle
import os
from Sat import SatManager as Sat
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins import V2
from postproc import masks
from commons.utils import addsep


def argument():
    parser = argparse.ArgumentParser(description = '''
    Apply check based on climatology to sat ORIG files
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--origdir', '-i',
                                type = str,
                                required = True,
                                help = ''' ORIG sat directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/ORIG/'''

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

    parser.add_argument(   '--checkdir', '-o',
                                type = str,
                                required = True,
                                help = ''' base of CHECKED and REJECTED directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/'''

                                )

    parser.add_argument(   '--statsdir', '-w',
                                type = str,
                                required = True,
                                help = ''' STATS  directory'''

                                )

    return parser.parse_args()


args = argument()


ORIGDIR   = addsep(args.origdir)
SUBMASKDIR  = addsep(args.submaskdir)
CHECKDIR  = addsep(args.checkdir)
STATSDIR  = addsep(args.statsdir)
CLIM_FILE = args.climfile

maskSat = getattr(masks,args.mesh)

reset = True #False

Timestart="19501231"
Time__end="20500101"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TL_orig = TimeList.fromfilenames(TI, ORIGDIR ,"*.nc",prefix='',dateformat='%Y%m%d')

somecheck = False
for filename in TL_orig.filelist:
    outfile = CHECKDIR + '/CHECKED/' + os.path.basename(filename)
    if not os.path.exists(outfile) :
        somecheck = True
        break

if somecheck or reset:
    print('Read climatology')
    MEAN,STD = Sat.readClimatology(CLIM_FILE)

    filemaskmed = SUBMASKDIR + 'maskmed_1km.npy'
    maskmed_1km = np.load(filemaskmed)
    numP = maskmed_1km.sum(axis=None)

    masksub_M = {}
    for sub in V2.P:
        filemasksub = SUBMASKDIR + 'masksub.' + sub.name + '.npy'
        masksub = np.load(filemasksub)
        masksub_M[sub.name] = np.zeros((maskSat.jpj,maskSat.jpi),dtype=bool)
        masksub_M[sub.name][maskmed_1km] = masksub


else:
    print('All checks done')



for iTime, filename in enumerate(TL_orig.filelist):
    outfile = CHECKDIR + '/CHECKED/' + os.path.basename(filename)
    exit_condition = os.path.exists(outfile) and (not reset)
    if exit_condition:
        continue
    iDate = TL_orig.Timelist[iTime] 
    date8 = '%04d' %iDate.year + '%02d' %iDate.month + '%02d' %iDate.day
    julian = int( iDate.strftime("%j") )
    print(' ... day ' + np.str(julian) + '  of ' + np.str(iDate.year))
    if julian == 366:
        julian = 365

    rejfile = CHECKDIR + '/REJECTED/' + date8 + '_rejected.nc'

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
    
    print 'Done check with ', filename, '  (',iTime+1,' of ', len(TL_orig.filelist), ')'
    print 'Rejection:  after check', counter_elim, ' values'
    print 'rejected for NAN in Climatology', counter_refNAN, ' values'
    Sat.dumpGenericNativefile(outfile, CHL_OUT, "CHL",mesh=maskSat)
    Sat.dumpGenericNativefile(rejfile, maskreject, "RejInd",mesh=maskSat)

    for sub in V2.P:
            masksubday = {}
            masksubday['ORIG'] = masksub_M[sub.name] & (cloudsLandTIME == False)
            masksubday['CHECK'] = masksub_M[sub.name] & (cloudsLandTIME == False) & (maskreject == 0)
            masksubday['ALLP'] = masksub_M[sub.name]

    for masktype in ['ORIG','CHECK','ALLP']:
        stats_clima = {}
        stats_day = {}
        print('    ---   Cycle on sub   ---')
        for sub in V2.P:
            stats_clima[sub.name] = np.zeros(12)
            stats_clima[sub.name][:] = np.nan
            stats_day[sub.name] = np.zeros(8)
            stats_day[sub.name][:] = np.nan

            climadmean = DAILY_REF_MEAN[masksubday[masktype]]
            climadmean[climadmean<0] = np.nan
            climadstd = DAILY_REF_STD[masksubday[masktype]]
            climadstd[climadmean<0] = np.nan

            climadmadd_2std = climadmean + 2*climadstd
            climadmsub_2std = climadmean - 2*climadstd

            chlsub = CHL_IN[masksubday[masktype]]
            chlsub[chlsub<0] = np.nan

            if np.nansum(masksubday[masktype])>0:
                #print(sub.name + ' ... Np ' + np.str(np.nansum(masksubday)))
                stats_clima[sub.name][0] = np.nanmean(climadmean)
                stats_clima[sub.name][1] = np.nanmin(climadmean)
                stats_clima[sub.name][2] = np.nanmax(climadmean)
                stats_clima[sub.name][3] = np.nanmean(climadmadd_2std)
                stats_clima[sub.name][4] = np.nanmin(climadmadd_2std)
                stats_clima[sub.name][5] = np.nanmax(climadmadd_2std)
                stats_clima[sub.name][6] = np.nanmean(climadmsub_2std)
                stats_clima[sub.name][7] = np.nanmin(climadmsub_2std)
                stats_clima[sub.name][8] = np.nanmax(climadmsub_2std)
                stats_clima[sub.name][9] = np.nanmean(climadstd)
                stats_clima[sub.name][10] = np.nanmin(climadstd)
                stats_clima[sub.name][11] = np.nanmax(climadstd)


                stats_day[sub.name][0] = np.nanmean(chlsub)
                stats_day[sub.name][1] = np.nanmin(chlsub)
                stats_day[sub.name][2] = np.nanmax(chlsub)
                stats_day[sub.name][3] = np.nanmean(np.abs(chlsub-climadmean)/climadstd)
                stats_day[sub.name][4] = np.nanmean(chlsub/climadmean)

                stats_day[sub.name][5] = np.nansum(masksubday['ORIG'])
                stats_day[sub.name][6] = np.nansum(outOfRange[masksubday['ORIG']])   #maskreject=1
                stats_day[sub.name][7] = np.nansum(~cloudsLandTIME & cloudlandsCLIM & masksubday['ORIG']) #maskreject=2


        fileclim = STATSDIR + 'stats_clim' + date8 + masktype + '.pkl'
        fid = open(fileclim,'wb')
        pickle.dump(stats_clima,fid)
        fid.close()

        filechlsub = STATSDIR + 'stats_time' + date8 + masktype + '.pkl'
        fid = open(filechlsub,'wb')
        pickle.dump(stats_day,fid)
        fid.close()

    print 'Done statistics with ', filename, '  (',iTime,' of ', len(TL_orig.filelist), ')'
    print '   ---------------------------------------------------   '


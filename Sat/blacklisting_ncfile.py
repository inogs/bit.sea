import numpy as np
import argparse
import pickle
import os
from Sat import SatManager as Sat
from Sat import interp2d
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons import timerequestors
from commons import netcdf3
from commons.mask import Mask
from basins import V2
from postproc import masks
from commons.utils import addsep
from commons.timeseries import TimeSeries


def argument():
    parser = argparse.ArgumentParser(description = '''
    Compose the blacklistinng nc file.
    Cycle on DA dates from archive.
    Dates of daily sat rejected (checked with avesat dates)
    Apply DA flag on all the dates to 1km resolution (interp2).
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--rejdir', '-r',
                                type = str,
                                required = True,
                                help = ''' Dir with .nc rejected from clima check'''

                                )

    parser.add_argument(   '--archivedir', '-d',
                                type = str,
                                required = True,
                                help = ''' Dir with .nc DA limitation'''
                                )

    parser.add_argument(   '--avesatdir', '-s',
                                type = str,
                                required = True,
                                help = ''' Dir with sat dates'''

                                )

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' Out dir with .nc file of rejected/DAlimited flags'''

                                )

    parser.add_argument(   '--satmesh', '-t',
                                type = str,
                                required = True,
                                choices = ['SatOrigMesh','V4mesh','V1mesh','KD490mesh','SAT1km_mesh', 'Mesh24'],
                                help = ''' Name of the mesh of sat ORIG'''
                                )

    parser.add_argument(   '--modmesh', '-m',
                                type = str,
                                required = True,
                                help = ''' Name of the mesh of model used in DA'''
                                )


    return parser.parse_args()


args = argument()


REJECTDIR  = addsep(args.rejdir)
ARCHIVEDIR = addsep(args.archivedir)
AVESATDIR  = addsep(args.avesatdir)
OUTDIR     = addsep(args.outdir)

maskSat = getattr(masks,args.satmesh)
maskMod = Mask(args.modmesh)

reset = True #False

Timestart="19501231"
Time__end="20170909"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")


TS = TimeSeries(TI,archive_dir=ARCHIVEDIR, \
                postfix_dir='DA__FREQ_1', \
                glob_pattern='limcorr.')
path_dadir = TS.get_runs([2]) # Tuesday assimilation
# TL_DA = TimeList.fromfilenames(TI, DADIR ,"*.txt", \
#                 prefix='',dateformat='%Y%m%d')



somecheck = False
for iDate,directory in path_dadir:
    date8DA = iDate.strftime('%Y%m%d')

    req = timerequestors.Weekly_req(iDate.year,iDate.month,iDate.day)
    TL_rej = TimeList.fromfilenames(req.time_interval, REJECTDIR,"*.nc", \
                prefix='',dateformat='%Y%m%d')
    for dd in TL_rej.Timelist:
        date8 = dd.strftime('%Y%m%d')
        outfile = OUTDIR + date8 + 'satflags.nc'
        if not os.path.exists(outfile):
            somecheck = True
            break
    if somecheck:
        break


if somecheck or reset:
    print('Some blacklisting file to be created')
    x = maskMod.xlevels[0,:]
    y = maskMod.ylevels[:,0]

    x1km = maskSat.lon
    y1km = maskSat.lat

    I_START, I_END = interp2d.array_of_indices_for_slicing(x, x1km)
    J_START, J_END = interp2d.array_of_indices_for_slicing(y, y1km)

    

else:
    print('All blacklisting files already done')



for iDate,directory in path_dadir:
    date8DA = iDate.strftime('%Y%m%d')
    fileDA = directory + '/limcorr.' + date8DA + '-000000.nc'
    
    file_weekdates = AVESATDIR + \
                     date8DA + 'weekdates.txt'
    weekdates = np.loadtxt(file_weekdates,dtype=str)
    Ndates = len(weekdates)

    req = timerequestors.Weekly_req(iDate.year,iDate.month,iDate.day)
    TL_rej = TimeList.fromfilenames(req.time_interval, REJECTDIR,"*_rejected.nc", \
                prefix='',dateformat='%Y%m%d')
    
    if (os.path.exists(fileDA)):
        flagDA = netcdf3.read_2d_file(fileDA,'lim_to_corr_FLAG')
        indlonflag,indlatflag = np.nonzero(flagDA)
        NflagDA = indlonflag.shape[0]
        maskDA_satmesh = np.zeros((maskSat.jpj,maskSat.jpi),dtype=int)
        for iiDA in range(NflagDA):
            maskDA_satmesh[I_START[indlonflag[iiDA]]:I_END[indlonflag[iiDA]],
                           J_START[indlatflag[iiDA]]:J_END[indlatflag[iiDA]]] = 1

    else:
        print(date8DA + ' Not existing flagDA file: ' + fileDA)
        continue
    
    print(iDate)


    flagmap = np.zeros((maskSat.jpj,maskSat.jpi))
    for iid,dd in enumerate(TL_rej.Timelist):
        date8 = dd.strftime('%Y%m%d')
        outfile = OUTDIR + date8 + 'satflags.nc'
        exit_condition = os.path.exists(outfile) and (not reset)
        if exit_condition:
            continue
        
        if Ndates==0:
            print(date8DA + ' file not used for avesat because Ndates<3')
            flagmap[:,:] = np.nan
        else:        
            if not(date8 in weekdates):
                print(date8 + ' in rejected but NOT in aveSat list: CHECK')
            fileflagclima = TL_rej.filelist[iid]
            flagclima = netcdf3.read_2d_file(fileflagclima,'RejInd')
            flagmap = flagclima
            flagmap[maskDA_satmesh==1] = 3

        Sat.dumpGenericNativefile(outfile, flagmap, "flagValues",mesh=maskSat)
    





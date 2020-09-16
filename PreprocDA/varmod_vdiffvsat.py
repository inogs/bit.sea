import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Calculate var_mod as difference between vardiff and varsat
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Outputdir'''
                            )

    parser.add_argument(   '--insat', '-s',
                            type = str,
                            required = True,
                            help = 'Input var sat'
                            )

    parser.add_argument(   '--indif', '-d',
                            type = str,
                            required = True,
                            help = 'Input var dif'
                            )

#    parser.add_argument(   '--limgib', '-g',
#                            type = str,
#                            required = True,
#                            help = 'Liongitude limit of Gibraltar'
#                            )

    parser.add_argument(   '--maskfile', '-m',
                            type = str,
                            required = True,
                            help = 'model maskfile meshmask.nc'
                            )

    parser.add_argument(   '--minvmod', '-v',
                            type = str,
                            required = True,
                            nargs=2,
                            help = 'min mod var quotas of vdiff in open sea and coast'
                            )

    return parser.parse_args()


args = argument()


import numpy as np
from commons.dataextractor import DataExtractor
from commons.mask import Mask 
from commons.utils import addsep
from postproc import masks
from Sat import SatManager

maskSat = getattr(masks,'Mesh24')

INSAT = addsep(args.insat)
INDIF = addsep(args.indif)
OUTDIR = addsep(args.outdir)
TheMask = Mask(args.maskfile)
minvmodLIST = args.minvmod

indGib,_ = TheMask.convert_lon_lat_to_indices(-5.2,35.9)

masksurf = TheMask.mask[0,:,:].copy()
masksurf[:,:indGib] = False

maskOpen = (TheMask.mask_at_level(200)) & (masksurf)
maskCoast =(TheMask.mask_at_level(200)==False) & (masksurf) 

minOpen = float(minvmodLIST[0])
minCoast = float(minvmodLIST[1])

for im in range(1,13):
    print 'Month %i' %im
    filesat = INSAT + "/var2Dsat.CCI.F7.2.%02d.nc"%(im)
    De = DataExtractor(TheMask,filename=filesat,varname='variance',dimvar=2)
    varSat = De.filled_values[:,:]
    filedif = INDIF + "/varErr.%02d.nc"%(im)
    De = DataExtractor(TheMask,filename=filedif,varname='variance',dimvar=2)
    varDif = De.filled_values[:,:]

    varMod = varDif-varSat

    maskminOpen = (varMod<minOpen*varDif) & maskOpen
    varMod[maskminOpen] = minOpen*varDif[maskminOpen]
    maskminCoast = (varMod<minCoast*varDif) & maskCoast
    varMod[maskminCoast] = minCoast*varDif[maskminCoast]

    iterfill = np.any(np.isnan(varMod[masksurf]))
    print '  N nan varMod points %i' %np.sum(np.isnan(varMod[masksurf]))
    if iterfill:
        for ifill in range(5):
            indxsp = np.argwhere(np.isnan(varMod) & masksurf)
            for jjp,iip in indxsp:
                varMod[jjp,iip] = np.nanmean(varMod[jjp-1:jjp+2,iip-1:iip+2])
            iterfnew = np.any(np.isnan(varMod[masksurf]))
            if iterfnew==False: 
                print "  Any nan after %i iteration" %(ifill+1)
                break
    if iterfnew:
        print '---- Still nan in variance with 5 iteration ----'

    MonthStr = '%02d'%(im)

    ratio = varMod/varSat
    fname = OUTDIR + '/ratiovarM_S.' + MonthStr + '.nc'
    print "\tsaving ", fname
    SatManager.dumpGenericNativefile(fname,ratio,'variance',maskSat)

    fname = OUTDIR + '/var2D.' + MonthStr + '.nc'
    print "\tsaving ", fname
    SatManager.dumpGenericNativefile(fname,varMod,'variance',maskSat)

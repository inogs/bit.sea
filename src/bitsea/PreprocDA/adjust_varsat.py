import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Adjust satellite variance files to have:
    i) larger variance in summer
    ii) smaller variance in coastal areas
    The changes can be done indipendently setting arguments to nan
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Outputdir'''
                            )

    parser.add_argument(   '--indir', '-i',
                            type = str,
                            required = True,
                            help = 'Input old VAR_SAT'
                            )

    parser.add_argument(   '--maskfile', '-m',
                            type = str,
                            required = True,
                            help = 'meshmask.nc file'
                            )

    parser.add_argument(   '--summerstart', '-s',
                            type = str,
                            required = True,
                            help = 'summer first month'
                            )

    parser.add_argument(   '--summerend', '-e',
                            type = str,
                            required = True,
                            help = 'summer last month'
                            )

    parser.add_argument(   '--additive', '-v',
                            type = str,
                            required = True,
                            help = 'additive error for summer (if nan: summer var maps unchanged)'
                            )

    parser.add_argument(   '--limitcoast', '-c',
                            type = str,
                            required = True,
                            help = 'Limit on variance in the coastal area (if nan: coastal areas unchanged)'
                            )

    parser.add_argument(   '--mesh', '-mm',
                            type = str,
                            required = True,
                            help = 'Mesh of sat map resolution'
                            )

    return parser.parse_args()

args = argument()
import numpy as np
# from bitsea.commons.Timelist import TimeList
# from bitsea.commons.time_interval import TimeInterval
# from bitsea.commons.timerequestors import Clim_month
from bitsea.commons.utils import addsep
from bitsea.commons.mask import Mask
from bitsea.postproc import masks
import bitsea.Sat.SatManager as Sat
import netCDF4



limerrcoast = args.limitcoast
if limerrcoast=='nan':
    limitcoast = False
else:
    limerrcoast = float(limerrcoast)
    limitcoast = True

summerstart = int(args.summerstart)
summerend = int(args.summerend)

additive = args.additive
if 'nan' in additive:
    adjustsummer = False
else:
    additive = float(additive)
    adjustsummer = True

MyMesh = getattr(masks,args.mesh)

TheMask = Mask.from_file(args.maskfile)

INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)


import os


mask200 = TheMask.mask_at_level(200)


def modify_summer(varsatin,additive):
    varsatout = (varsatin**.5 + additive)**2
    return varsatout

def modify_coast(varsatin,limerrcoast):
    masklim = varsatin>limerrcoast**2
    varsatout = varsatin.copy()
    varsatout[masklim & mask200==False] = limerrcoast**2
    return varsatout

def saveNC(variance,filevarnew):
    print(f' ---- Saving {filevarnew} ----')
    VV = netCDF4.Dataset(filevarnew,'w')
    VV.createDimension('lon', MyMesh.jpi)
    VV.createDimension('lat', MyMesh.jpj)
    ncvar = VV.createVariable('variance', 'f', ('lat', 'lon'))
    setattr(ncvar,'missing_value',1.e+20)
    ncvar[:] = variance
    ncvar = VV.createVariable('lat','f',('lat',))                          
    ncvar[:] = MyMesh.lat                                                       
    ncvar = VV.createVariable('lon','f',('lon',))                          
    ncvar[:] = MyMesh.lon                                                       
    VV.close()     
    return



for mm in range(1,13):
    filevarname = f'var2D.{mm:02d}.nc'
    print(f'---------- {filevarname} ----------')
    if (limitcoast==False) & (adjustsummer==False):
        print('not modifying variance')
        os.system(f'cp {INDIR}/{filevarname} {OUTDIR}/{filevarname}')
        continue
    if limitcoast:
        print('modify coastal variance')
        VV = netCDF4.Dataset(f'{INDIR}/{filevarname}','r')
        varold = VV.variables['variance'][:].data
        VV.close()
        varnew = modify_coast(varold,limerrcoast)
        if (mm>=summerstart and mm<=summerend) & (adjustsummer):
            print('modify summer variance')
            varnew = modify_summer(varnew,additive)
        fileout = f'{OUTDIR}/{filevarname}'
        saveNC(varnew,f'{OUTDIR}/{filevarname}')
    else:
        if adjustsummer:
            if (mm>=summerstart and mm<=summerend):
                VV = netCDF4.Dataset(f'{INDIR}/{filevarname}','r')
                varold = VV.variables['variance'][:].data
                VV.close()
                print('modify summer variance')
                varnew = modify_summer(varold,additive)
                fileout = f'{OUTDIR}/{filevarname}'
                saveNC(varnew,f'{OUTDIR}/{filevarname}')
            elif (mm<summerstart or mm>summerend):
                print('not modifying summer variance')
                os.system(f'cp {INDIR}/{filevarname} {OUTDIR}/{filevarname}')    



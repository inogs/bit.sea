import argparse
from numpy import dtype
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates png maps, based on matplotlib.imshow(), of time averaged fields
    
    Brief description of the algorithm:
    INPUTS                                               OPERATOR              OUTPUT
    (InputDir, varname, StartTime, EndTime)  ---->     time average       ---> 3D field
    (3d_field, Layer list )                  ----> vertical mean/integral ---> 2d map
    
    
    Example of output file: 
    Map_pCO2_Ave.2016-2015_Ave0060-0100m-0060m.png
    
    Caveats:
    A unit conversion is performed, about ppn.
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Output dir'''
                            )

    parser.add_argument(   '--inputdir', '-i',
                            type = str,
                            required = True,
                            default = './',
                            help = 'Input dir')

    parser.add_argument(   '--varname', '-v',
                                type = str,
                                required = True,
                                default = '',
                                choices = ['P_l','P_i','N1p', 'N3n', 'O2o', 'pCO2','PH','pH','ppn','P_c','Ac','DIC','Z_c'] )
    parser.add_argument(   '--plotlistfile', '-l',
                                type = str,
                                required = True,
                                help = '''Plotlist_bio.xml"
                                ''')
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = 'Path of the mask file')
    parser.add_argument(   '--optype', '-t',
                                type = str,
                                required = True,
                                default = '',
                                choices = ['integral','mean'],
                                help ="  INTEGRALE:  * heigth of the layer, MEDIA    :  average of layer")
    parser.add_argument(   '--starttime','-s',
                                type = str,
                                required = True,
                                help = 'start date in yyyymmdd format')
    parser.add_argument(   '--endtime','-e',
                                type = str,
                                required = True,
                                help = 'start date in yyyymmdd format')      

    return parser.parse_args()
args = argument()

import numpy as np

from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.commons.mask import Mask
from bitsea.commons.layer import Layer
from bitsea.layer_integral.mapbuilder import MapBuilder, Plot
from bitsea.layer_integral.mapplot import mapplot,pl
from bitsea.commons.dataextractor import DataExtractor
from bitsea.commons.time_averagers import TimeAverager3D
from bitsea.layer_integral import coastline
import bitsea.commons.timerequestors as requestors
from bitsea.commons.utils import addsep
from bitsea.commons.xml_module import *
from xml.dom import minidom
from bitsea.commons import netcdf3

xmldoc = minidom.parse(args.plotlistfile)


INPUTDIR  = addsep(args.inputdir)
OUTPUTDIR = addsep(args.outdir)
var       = args.varname


for lm in xmldoc.getElementsByTagName("LayersMaps"):
    for pdef in get_subelements(lm, "plots"):
        filevar = get_node_attr(pdef, "var")        
        if not filevar == var : continue
        PLOT = Plot(get_node_attr(pdef, "var"), get_node_attr(pdef, "longname"), get_node_attr(pdef, "plotunits"), [], [0,1] )
        for d in get_subelements(pdef, "depth"):
            levplot  = float(get_node_attr(d, "levplot")) 
            clim     = eval(get_node_attr(d, "clim"))
            L = Layer(get_node_attr(d,"top"), get_node_attr(d, "bottom"))
            PLOT.append_layer(L,clim=clim, mapdepthfilter=levplot)
            




clon,clat = coastline.get()
TheMask = Mask.from_file(args.maskfile)

CONVERSION_DICT={
         'ppn' : 365./1000,
         'Z_c' : 1./1000,
         'O2o' : 1,
         'N1p' : 1,
         'N3n' : 1,
         'PH'  : 1,
         'pH'  : 1, 
         'pCO2': 1,
         'P_l' : 1,
         'P_c' : 1., # 1./1000, #divided by 1000 in the integral
         'P_i' : 1,
	     'Ac'  : 1,
         'DIC' : 1
         }

MONTH_STRING = ["January","February","March","April","May","June","July","August","September","October","November","December"]
TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")
req_label = "Ave." + str(TI.start_time.year) + "-" +str(TI.end_time.year-1)

TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*.nc",filtervar=var)
if TL.inputFrequency is None:
    TL.inputFrequency='monthly'
    print ("inputFrequency forced to monthly because of selection of single time")

req = requestors.Generic_req(TI)
indexes,weights = TL.select(req)


VARCONV=CONVERSION_DICT[var]
# setting up filelist for requested period -----------------
filelist=[]
for k in indexes:
    t = TL.Timelist[k]
    filename = INPUTDIR + "ave." + t.strftime("%Y%m%d-%H:%M:%S") + "." + var + ".nc"
    filelist.append(filename)
# ----------------------------------------------------------
print ("time averaging ...")
if (var == "CO2airflux"):
   M3d     = TimeAverager2D(filelist, weights, var, TheMask)
else:
   M3d     = TimeAverager3D(filelist, weights, var, TheMask)
print ("... done.")
for il, layer in enumerate(PLOT.layerlist):
    z_mask = PLOT.depthfilters[il]
    z_mask_string = "-%04gm" %z_mask
    if  args.optype=='integral':
        outfile    = OUTPUTDIR + "Map_" + var + "_" + req_label + "_Int" + layer.longname() + z_mask_string  + ".png"
    else:
        outfile    = OUTPUTDIR + "Map_" + var + "_" + req_label + "_Ave" + layer.longname() + z_mask_string  + ".png"

    De      = DataExtractor(TheMask,rawdata=M3d)

    if args.optype=='integral':
        integrated = MapBuilder.get_layer_integral(De, layer)
    else:
        integrated = MapBuilder.get_layer_average(De, layer)
    integrated=integrated * VARCONV


    mask=TheMask.mask_at_level(z_mask)
#        clim = [M3d[TheMask.mask].min(), M3d[TheMask.mask].max()]
    clim=PLOT.climlist[il]
    integrated_masked=integrated*mask # taglio il costiero
    integrated_masked[integrated_masked==0]=np.nan # sostituisco gli 0 con i NAN

    fig,ax     = mapplot({'clim':clim, 'data':integrated_masked, },fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
    ax.set_xlim([-5,36])
    ax.set_ylim([30,46])
    ax.set_xlabel('Lon').set_fontsize(11)
    ax.set_ylabel('Lat').set_fontsize(12)
#    ax.ticklabel_format(fontsize=10)
    ax.tick_params(axis='x', labelsize=10)
#    ax.text(-4,44.5,PLOT.name() + ' [' + PLOT.units() + ']',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')
    if (var=="Z_c"):
       ax.text(-4,44.5, 'ZOOPL [' + PLOT.units() + ']',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')
    else:
       ax.text(-4,44.5,var + ' [' + PLOT.units() + ']',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')

    ax.xaxis.set_ticks(np.arange(-2,36,6))
    ax.yaxis.set_ticks(np.arange(30,46,4))
    #ax.text(-4,30.5,req_label,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
    ax.grid()
    if (var=="Z_c"):
       title = "%s %s %s" % ('annual', 'ZOOPL', layer.__repr__())
    else:
        title = "%s %s %s" % ('annual', var, layer.__repr__())
#    if (var == "PH"): title = "%s %s %s" % ('annual', "pH$\mathrm{_T}$", layer.__repr__())
#    title = "%s %s %s" % (MONTH_STRING[TI.start_time.month - 1], var, layer.__repr__())
    if (var == "pH"): title = "%s %s %s" % (MONTH_STRING[TI.start_time.month - 1], "pH$\mathrm{_T}$", layer.__repr__())
    if (var == "pCO2"): title = "%s %s %s" % (MONTH_STRING[TI.start_time.month - 1], "pCO2", layer.__repr__())
    if (var=="P_c"): title = "%s %s %s" % ('annual', "phytoplankton biomass", layer.__repr__())
    fig.suptitle(title)
    fig.savefig(outfile)
    pl.close(fig)
    if (var == "ppn"): 
        ncfile = OUTPUTDIR + "Map_" + var + "_" + req_label + "_Int" + layer.longname() + z_mask_string  + ".nc"
#    netcdf3.write_2d_file(integrated_masked,"ppn",ncfile,mask)
        netcdf3.write_2d_file(integrated_masked,"ppn",ncfile,TheMask)
        
    if (var == "ppn"or var == "Z_c" or var == "P_c"):
        from bitsea.basins import V2 as OGS
        from bitsea.commons.submask import SubMask
        from bitsea.basins.basin import ComposedBasin
        from bitsea.commons.utils import writetable
        tablefile = OUTPUTDIR + '/' + var + '_mean_basin.txt'
        OGSred = ComposedBasin('OGSred',[OGS.alb, \
             #       OGS.swm1, OGS.swm2, OGS.nwm, OGS.tyr, \
                    OGS.swm, OGS.nwm, OGS.tyr, \
                    OGS.adr, OGS.aeg, OGS.ion, OGS.lev , OGS.med], \
                   # 'Gruped Subbasin for ppn analysis')
                    'Gruped Subbasin for integral variable analysis')
        SUBlist = OGSred.basin_list
        nSub   = len(OGSred.basin_list)
        rows_names=[sub.name for sub in SUBlist]
        ppn_submean = np.zeros((nSub,2),np.float32)*np.nan # FIRST COL: mean; SECOND COL: std
        for isub, sub in enumerate(OGSred):
            S = SubMask(sub, TheMask)
            mask2d=S.mask[0,:,:]
            ppn_submean[isub,0] = integrated_masked[mask2d].mean()
            ppn_submean[isub,1] = integrated_masked[mask2d].std()

#        writetable(tablefile,ppn_submean,rows_names,['mean'])
        writetable(tablefile,ppn_submean,rows_names,['MEAN' ,'STD'],fmt=("%8.2f\t %8.2f"))

        

#        ncfile = OUTPUTDIR + "Map_" + var + "_" + req_label + "_Int" + layer.longname() + z_mask_string  + ".nc"
#        netcdf3.write_2d_file(integrated_masked,"Z_c",ncfile,TheMask)
#
#        from bitsea.basins import V2 as OGS
#        from bitsea.commons.submask import SubMask
#        from bitsea.basins.basin import ComposedBasin
##        from bitsea.commons.utils import writetable
#        OGSred = ComposedBasin('OGSred',[OGS.alb, OGS.swm1, OGS.swm2, OGS.nwm, OGS.tyr, OGS.adr, OGS.ion, OGS.lev , OGS.med],'Gruped Subbasin for ppn analysis')
#        SUBlist = OGSred.basin_list
#        nSub   = len(OGSred.basin_list)
#        rows_names=[sub.name for sub in SUBlist]
#column_names=['bias','rmse','corr']
#        ppn_submean = np.zeros((nSub,1),np.float32)*np.nan
#        for isub, sub in enumerate(OGSred):
#            S = SubMask(sub, TheMask)
#            mask2d=S.mask[0,:,:]
#            ppn_submean[isub,0] = integrated_masked[mask2d].mean()
#
#        writetable(OUTPUTDIR + var + "MEAN-BASIN.txt", ppn_submean, rows_names, ['mean'])
#  

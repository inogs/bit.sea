import argparse
from bitsea.utilities.argparse_types import date_from_str
from bitsea.utilities.argparse_types import existing_dir_path, existing_file_path
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates png maps, based on matplotlib.imshow(), of time averaged fields
    
    Brief description of the algorithm:
    INPUTS                                               OPERATOR              OUTPUT
    (InputDir, varname, StartTime, EndTime)  ---->     time average       ---> 3D field
    (3d_field, Layer list )                  ----> vertical mean/integral ---> 2d map
    
    
    Example of output file: 
    pCO2_2016-2015_Ave0060-0100m-0060m.png
    
    Caveats:
    A unit conversion is performed, about ppn.
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = existing_dir_path,
                            required =True,
                            help = ''' Output dir'''
                            )

    parser.add_argument(   '--inputdir', '-i',
                            type = existing_dir_path,
                            required = True,
                            help = 'Input dir')

    parser.add_argument(   '--varname', '-v',
                                type = str,
                                required = True,
                                default = '',
                                choices = ['P_l','Z_c','N1p', 'N3n', 'O2o', 'pCO2','PH','pH','ppn','P_c','Ac','ALK','DIC','netPPYc','N4n','N5s','CO2airflux'] )
    parser.add_argument(   '--plotlistfile', '-l',
                                type = str,
                                required = True,
                                help = '''Plotlist_bio.xml"
                                ''')
    parser.add_argument(   '--maskfile', '-m',
                                type = existing_file_path,
                                required = True,
                                help = 'Path of the mask file')
    parser.add_argument(   '--optype', '-t',
                                type = str,
                                required = True,
                                default = '',
                                choices = ['integral','mean'],
                                help ="  INTEGRALE:  * heigth of the layer, MEDIA    :  average of layer")
    parser.add_argument(   '--starttime','-s',
                                type = date_from_str,
                                required = True,
                                help = 'start date in yyyymmdd format')
    parser.add_argument(   '--endtime','-e',
                                type = date_from_str,
                                required = True,
                                help = 'start date in yyyymmdd format')      

    return parser.parse_args()
args = argument()

import numpy as np
import matplotlib
matplotlib.use('Agg')
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.Timelist import TimeList
from bitsea.commons.mask import Mask
from bitsea.commons.layer import Layer
from bitsea.layer_integral.mapbuilder import MapBuilder, Plot
from bitsea.layer_integral.mapplot import mapplot,pl
from bitsea.layer_integral.mapplot import mapplotlog
from bitsea.commons.dataextractor import DataExtractor
from bitsea.commons.time_averagers import TimeAverager3D
from bitsea.layer_integral import coastline
import bitsea.commons.timerequestors as requestors
from bitsea.commons.xml_module import get_subelements, get_node_attr
from xml.dom import minidom
from bitsea.commons import netcdf3

xmldoc = minidom.parse(args.plotlistfile)


INPUTDIR  = args.inputdir
OUTPUTDIR = args.outdir
var       = args.varname

if  args.optype=='integral':
    maptype="Int"
else:
    maptype="Ave"

varlog = ''
if var == 'P_l': varlog = 'P_llog'

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
            
for lm in xmldoc.getElementsByTagName("LayersMaps"):
    for pdef in get_subelements(lm, "plots"):
        filevar = get_node_attr(pdef, "var")        
        if not filevar == varlog : continue
        PLOTlog = Plot(get_node_attr(pdef, "var"), get_node_attr(pdef, "longname"), get_node_attr(pdef, "plotunits"), [], [0,1] )
        for d in get_subelements(pdef, "depth"):
            levplot  = float(get_node_attr(d, "levplot")) 
            clim     = eval(get_node_attr(d, "clim"))
            L = Layer(get_node_attr(d,"top"), get_node_attr(d, "bottom"))
            PLOTlog.append_layer(L,clim=clim, mapdepthfilter=levplot)



clon,clat = coastline.get()
TheMask = Mask.from_file(args.maskfile)

CONVERSION_DICT={
     'netPPYc' : 365./1000,
         'ppn' : 365./1000,
         'O2o' : 1,
         'N1p' : 1,
         'N3n' : 1,
         'N4n' : 1,
         'N5s' : 1,
         'PH'  : 1,
         'pH'  : 1, 
         'pCO2': 1,
         'P_l' : 1,
         'P_c' : 1,
         'Z_c' : 1,
	 'Ac'  : 1,
         'ALK' : 1,
         'DIC' : 1, 
  'CO2airflux' : 1
         }

MONTH_STRING = ["January","February","March","April","May","June","July","August","September","October","November","December"]
TI = TimeInterval.fromdatetimes(args.starttime,args.endtime)
if args.endtime.year > args.starttime.year:
    req_label=f"{args.starttime.year}_{args.endtime.year}"
else:
    req_label=f"{args.starttime.year}"

TL = TimeList.fromfilenames(TI, INPUTDIR,"ave*.nc",filtervar="."+var)
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
    filename = INPUTDIR / ("ave." + t.strftime("%Y%m%d-%H:%M:%S") + "." + var + ".nc")
    filelist.append(filename)
# ----------------------------------------------------------
print ("time averaging ...")
M3d     = TimeAverager3D(filelist, weights, var, TheMask)
print ("... done.")
for il, layer in enumerate(PLOT.layerlist):
    z_mask = PLOT.depthfilters[il]
    z_mask_str = "-%04gm" %z_mask
    Lstr=layer.longname()

    filename = f"{var}_{req_label}_{maptype}{Lstr}{z_mask_str}"
    outfile=OUTPUTDIR / f"{filename}.png"
    outfilelog= OUTPUTDIR / f"Log{filename}.png"


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

    fig,ax     = mapplot({'clim':clim, 'data':integrated_masked, }, \
        fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
    ax.set_xlim([-6,36])
    ax.set_ylim([30,46])
    ax.set_xlabel('Lon').set_fontsize(11)
    ax.set_ylabel('Lat').set_fontsize(12)
    #ax.ticklabel_format(fontsize=10)
    ax.tick_params(axis='x', labelsize=10)
    ax.text(-4,44.5,var + ' [' + PLOT.units() + ']',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')

    ax.xaxis.set_ticks(np.arange(-6,36,6))
    ax.yaxis.set_ticks(np.arange(30,46,4))
    #ax.text(-4,30.5,req_label,horizontalalignment='left',verticalalignment='center',fontsize=13, color='black')
    ax.grid()
    title = "%s %s %s" % (req_label, var, layer.__repr__())
    fig.suptitle(title)
    fig.savefig(outfile)
    pl.close(fig)
    if (var == "ppn"):
        ncfile = OUTPUTDIR / f"{filename}.nc"
        netcdf3.write_2d_file(integrated_masked,"ppn",ncfile,TheMask)

    if (var == 'P_l'):
        climlog=PLOTlog.climlist[il]
        fig,ax = mapplotlog({'clim':climlog, 'data':integrated_masked, }, \
            fig=None,ax=None,mask=TheMask,coastline_lon=clon,coastline_lat=clat)
        ax.set_xlim([-6,36])
        ax.set_ylim([30,46])
        ax.set_xlabel('Lon').set_fontsize(11)
        ax.set_ylabel('Lat').set_fontsize(12)
        # ax.ticklabel_format(fontsize=10)
        ax.tick_params(axis='x', labelsize=10)
        ax.text(-4,44.5,var + ' [' + PLOT.units() + ']',horizontalalignment='left',verticalalignment='center',fontsize=14, color='black')

        ax.xaxis.set_ticks(np.arange(-6,36,6))
        ax.yaxis.set_ticks(np.arange(30,46,4))
        ax.grid()
        title = "%s %s %s" % (req_label, var, layer.__repr__())
        fig.suptitle(title)
        fig.savefig(outfilelog)
        pl.close(fig)

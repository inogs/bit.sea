import matplotlib.pyplot as pl
import numpy as np
import matplotlib.patches as mpatches
import scipy.io.netcdf as NC
#import libxmp, libxmp.utils
#from libxmp import XMPFiles, consts
from layer_integral import coastline

def figure_generator(p):
    ''' Generates a figure to plot the matchups related to a bioFloat cycle
    There are 6 axes: map, temperature, salinity, chl, oxygen and nitrate
    plus 4 new axes related to PAR, P_c, POC, pH

    Arguments:
    * p * is a profile object

    Returns
    fig, axes (array of axes handlers)
    '''
    fig, axs = pl.subplots(2,5, facecolor='w', edgecolor='k')
    hsize=16
    vsize=12
    fig.set_size_inches(hsize,vsize)
    #figtitle = " date="+p.time.strftime('%Y/%m/%d')+" float="+p.name()
    #fig.set_title(figtitle)
    fig.subplots_adjust(hspace = 0.3, wspace=0.3)
    axs = axs.ravel()

    ax = axs[0]
    c_lon, c_lat=coastline.get()
    ax.plot(c_lon,c_lat, color='#000000',linewidth=0.5)
    ax.plot(p.lon,p.lat,'ro')
    ax.set_xticks(np.arange(-6,36,2))
#    ax.set_yticks(np.arange(-30,46))
    ax.set_yticks(np.arange(0,100,2))
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    ax.set_title(p.time.strftime('%Y/%m/%d'))
    extent=10 #6 #degrees
    ax.set_xlim([p.lon -extent/2, p.lon+extent/2])
    ax.set_ylim([p.lat -extent/2, p.lat+extent/2])
    bbox=ax.get_position()

    deltax, _ =bbox.size
    new_deltay = deltax* hsize/vsize
    bottom = bbox.ymax - new_deltay
    ax.set_position([bbox.xmin, bottom, deltax, new_deltay])

    floatlabel = 'Float \n'+ p.name() +" - "+str(p._my_float.cycle)
    b_patch = mpatches.Patch(color='red', label='Model')
#    g_patch = mpatches.Patch(color='blue', label='Float')
    g_patch = mpatches.Patch(color='blue', label=floatlabel)
    ax.legend(handles=[b_patch,g_patch], bbox_to_anchor=(0, -0.5), loc=2)

    for ax in axs[1:]:
        ax.set_ylim(0,400)
        ax.locator_params(axis='x',nbins=4)
        ax.yaxis.grid()

    for ax in [axs[2], axs[3], axs[4], axs[6], axs[7], axs[8], axs[9]]:
        ax.set_yticklabels([])

    return fig,axs



def ncwriter(filenc,zlevels_out,profileobj):
    ''' Generates a NetCDF file of matchups related to a bioFloat cycle
    There are 6 axes: map, temperature, salinity, chl, oxygen and nitrate

    Arguments:
    * filenc      * is a output file name
    * zlevels_out * is the array of depths
    * profileobj  * is a profile object

    Returns
    f, model_handlers, float_handlers
    '''

    depths = len(zlevels_out)
    f = NC.netcdf_file(filenc, 'w')
    f.createDimension('levels', depths)
    f.createDimension('pos', 1)
    lon=f.createVariable('longitude', 'f', ('pos',))
    lat=f.createVariable('latitude' , 'f', ('pos',))
    lon[:] = profileobj.lon
    lat[:] = profileobj.lat
    m_array = f.createVariable('lev', 'f', ('levels',))
    m_array[:] = zlevels_out[:]
    setattr(f, 'time', profileobj.time.strftime("%Y%m%d-%H:%M:%S"))

    model_handlers=[]
    float_handlers=[]
    for model_varname in ['P_l','O2o','N3n','votemper','vosaline','EIR','POC',"P_c", "pH" ]:
        name_var = model_varname+"_model"
        m_array = f.createVariable(name_var, 'f', ('levels',))
        setattr(m_array, 'missing_value', 1.e+20)
        setattr(m_array, 'fillValue', 1.e+20)
        m_array[:] = np.ones((depths,),np.float32)*1.e+20

        name_var = model_varname+"_float"
        f_array = f.createVariable(name_var, 'f', ('levels',))
        setattr(f_array, 'missing_value', 1.e+20)
        setattr(f_array, 'fillValue', 1.e+20)
        f_array[:] = np.ones((depths,),np.float32)*1.e+20
        model_handlers.append(m_array)
        float_handlers.append(f_array)
    return f, model_handlers, float_handlers




#def add_metadata(filepng,p):
#    xmpfile = XMPFiles( file_path=filepng, open_forupdate=True )
#    xmp = xmpfile.get_xmp()
#    if xmp is None:
#        xmp = libxmp.XMPMeta()
#    xmp.set_property(consts.XMP_NS_DC, 'float', p.name() )
#    xmp.set_property(consts.XMP_NS_DC, 'date', p.time.strftime('%Y%m%d') )
#    xmp.set_property(consts.XMP_NS_DC, 'hour', p.time.strftime('%H:%M:%S') )
#    xmp.set_property(consts.XMP_NS_DC, 'position.lat',str(p.lat)+"N")
#    xmp.set_property(consts.XMP_NS_DC, 'position.lon',str(p.lon)+"E")
#    xmpfile.put_xmp(xmp)
#    xmpfile.close_file()

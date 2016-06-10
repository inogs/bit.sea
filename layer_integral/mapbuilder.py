# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import os
import warnings
import numpy as np
from xml.dom import minidom
from ast import literal_eval

from commons.mask import Mask
from commons.layer import Layer
from commons.utils import is_number, get_date_string
from commons.xml_module import *
from commons.dataextractor import DataExtractor
from commons.dataextractor import NotFoundError
from mapplot import mapplot,mapplot_medeaf,mapplot_nocolor
import datetime
import pylab as pl

def warn_user(msg):
    warnings.warn(msg, SyntaxWarning, stacklevel=2)

class Plot(object):
    def __init__(self, varname, longvarname, units, layerlist, clim):
        #Input validation
        self.__varname = str(varname)
        self.__longvarname = str(longvarname)
        print "units = ", units
        self.__units = str(units.encode('utf-8'))
        if not isinstance(layerlist, (list, tuple)) or ((len(layerlist) > 0) and not isinstance(layerlist[0], (Layer,))):
            raise ValueError("layerlist must be a list of Layers")
        self.__layerlist = layerlist
        if not (clim is None):
            if not isinstance(clim, (list, tuple)) or (len(clim) != 2) or not (is_number(clim[0]) and is_number(clim[1])):
                raise ValueError("clim must be a list of two numbers")
        self.__clim = clim
        self.__climslist = list()

    @property
    def varname(self):
        return self.__varname
    def longvarname(self):
        return self.__longvarname
    def units(self):
        return self.__units

    @property
    def layerlist(self):
        return self.__layerlist

    @property
    def clim(self):
        return self.__clim

    @property
    def climlist(self):
        return self.__climslist

    def append_layer(self, layer, clim=None):
        if not isinstance(layer, (Layer,)):
            raise ValueError("layer must be a Layer object")
        self.__layerlist.append(layer)
        self.__climslist.append(clim)

class MapBuilder(object):

    def __init__(self, plotlistfile, netcdffileslist, maskfile, outputdir):
        if isinstance(netcdffileslist, (list, tuple)):
            self.__netcdffileslist = netcdffileslist
        else:
            self.__netcdffileslist = [ str(netcdffileslist) ]
        if os.path.isdir(outputdir):
            self.__outputdir = outputdir
        else:
            raise ValueError("outputdir must be a path to a directory")
        self._mask = Mask(maskfile)
        xmldoc = minidom.parse(plotlistfile)
        self.__plotlist = list()
        #For each LayersMaps element
        for lm in xmldoc.getElementsByTagName("LayersMaps"):
            #For each plots element
            for pdef in get_subelements(lm, "plots"):
                clim = get_node_attr(pdef, "clim")
                if not (clim is None):
                    clim = literal_eval(clim)
                plot = Plot(get_node_attr(pdef, "var"), get_node_attr(pdef, "longname"), get_node_attr(pdef, "units"), [], clim)
                #For each depth element
                for d in get_subelements(pdef, "depth"):
                    clim = get_node_attr(d, "clim")
                    if not (clim is None):
                        clim = literal_eval(clim)
                    plot.append_layer(Layer(get_node_attr(d, "top"), get_node_attr(d, "bottom")), clim)
                self.__plotlist.append(plot)

    def plot_maps_data(self, coastline_lon=None, coastline_lat=None, maptype=0):
        '''
        Generator of a large set of images.
        Arguments : 
          * coastline_lon *  1-D array 
          * coastline_lat *  1-D array
          * map_type      *  integer, flag used to choose a plot method in mapplot module
                = 0, the default, to call mapplot.mapplot()
                = 1, to call mapplot.mapplot_onlycolor()
                = 2, to call mapplot.mapplot_nocolor()
        
        '''
        fig = None
        ax = None
        for f in self.__netcdffileslist:
            longdate , shortdate = get_date_string(f)
            for p in self.__plotlist:
                try:
                    de = DataExtractor(self._mask, filename=f, varname=p.varname)
                except NotFoundError as e:
                    msg="File: %s\n%s" % (f, e)
                    warn_user(msg)
                    continue
                for i,l in enumerate(p.layerlist):
                    outfile = "%s/ave.%s.%s.%s" % (self.__outputdir,shortdate, p.varname, l)
                    mapdata = MapBuilder.get_layer_average(de, l)
                    try:
                        clim = p.climlist[i]
                        if clim is None:
                            raise ValueError
                    except ValueError:
                        if p.clim is None:
                            raise ValueError("No clim defined for %s %s" % (p.varname, repr(l)))
                        else:
                            clim = p.clim
                    if maptype == 0:
                        fig, ax = mapplot({'varname':p.varname, 'clim':clim, 'layer':l, 'data':mapdata, 'date':longdate}, fig=fig, ax=ax, mask=self._mask, ncolors=24, coastline_lon=coastline_lon, coastline_lat=coastline_lat)
                        fig.savefig(outfile + ".png")
                    if maptype == 1:
                        dateobj=datetime.datetime.strptime(shortdate,'%Y%m%d')
                        mapdict={'varname':p.varname, 'longname':p.longvarname(), 'clim':clim, 'layer':l, 'data':mapdata, 'date':dateobj,'units':p.units()}
                        fig, ax = mapplot_medeaf(mapdict, fig=fig, ax=ax, mask=self._mask, ncolors=24)
                        fig.savefig(outfile + ".png",dpi=86)
                        pl.close(fig)
                    if maptype == 2:
                        fig, ax = mapplot_nocolor({'varname':p.varname, 'clim':clim, 'layer':l, 'data':mapdata, 'date':longdate}, fig=fig, ax=ax, mask=self._mask, ncolors=24, coastline_lon=coastline_lon, coastline_lat=coastline_lat)
                        fig.canvas.print_figure(outfile + ".svg")
                        return
                    
                fig = None
                ax = None

    @staticmethod
    def get_layer_average(data_extractor, layer):
        """Returns a 2D NumPy array with the average weighted over depth.

        Args:
            - *data_extractor*: a DataExtractor object that contains the data to
              extract.
            - *layer*: a Layer object that contains the vertical boundaries for
              the computation.
        """
        if not isinstance(layer, (Layer,)):
            raise ValueError("layer must be a Layer object")
        if not isinstance(data_extractor, (DataExtractor,)):
            raise ValueError("dataextractor must be a DataExtractor object")
        #Find Z indices
        top_index = np.where(data_extractor._mask.zlevels >= layer.top)[0][0]
        bottom_index = np.where(data_extractor._mask.zlevels < layer.bottom)[0][-1]
        if top_index == bottom_index:
            #Just one layer so we return the sliced data
            output = data_extractor.get_filled_values[top_index,:,:]
            return output
        #Workaround for Python ranges
        bottom_index += 1
        #Build local mask matrix
        lmask = np.array(data_extractor._mask.mask[top_index:bottom_index,:,:], dtype=np.double)
        #Build dz matrix
        dzm = np.ones_like(lmask, dtype=np.double)
        j = 0
        for i in range(top_index, bottom_index):
            dzm[j,:,:] = data_extractor._mask.dz[i]
            j += 1
        #Get the slice of the values
        v = np.array(data_extractor.values[top_index:bottom_index,:,:])
        #Build integral matrix (2D)
        integral = (v * dzm * lmask).sum(axis=0)
        #Build height matrix (2D)
        height = (dzm * lmask).sum(axis=0)
        indexmask = [height > 0]
        #Build output matrix (2D)
        output = np.full_like(integral, data_extractor.fill_value, dtype=np.double)
        output[indexmask] = integral[indexmask] / height[indexmask]
        return output


    @staticmethod
    def get_layer_integral(data_extractor,layer):
        """Returns a 2D NumPy array with the integral over depth.

        Args:
            - *data_extractor*: a DataExtractor object that contains the data to
              extract.
            - *layer*: a Layer object that contains the vertical boundaries for
              the computation.
        """
        if not isinstance(layer, (Layer,)):
            raise ValueError("layer must be a Layer object")
        if not isinstance(data_extractor, (DataExtractor,)):
            raise ValueError("dataextractor must be a DataExtractor object")
        #Find Z indices
        top_index = np.where(data_extractor._mask.zlevels >= layer.top)[0][0]
        bottom_index = np.where(data_extractor._mask.zlevels < layer.bottom)[0][-1]
        if top_index == bottom_index:
            #Just one layer so we return the sliced data
            output = data_extractor.get_filled_values[top_index,:,:]
            return output
        #Workaround for Python ranges
        bottom_index += 1
        #Build local mask matrix
        lmask = np.array(data_extractor._mask.mask[top_index:bottom_index,:,:], dtype=np.double)
        #Build dz matrix
        dzm = np.ones_like(lmask, dtype=np.double)
        j = 0
        for i in range(top_index, bottom_index):
            dzm[j,:,:] = data_extractor._mask.dz[i]
            j += 1
        #Get the slice of the values
        v = np.array(data_extractor.values[top_index:bottom_index,:,:])
        #Build integral matrix (2D)
        integral = (v * dzm * lmask).sum(axis=0)
        #Build height matrix (2D)
        height = (dzm * lmask).sum(axis=0)
        indexmask = [height > 0]
        #Build output matrix (2D)
        output = np.full_like(integral, data_extractor.fill_value, dtype=np.double)
        output[indexmask] = integral[indexmask]
        return output
        
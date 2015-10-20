# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import os
import re
import numpy as np
from xml.dom import minidom
from ast import literal_eval

from commons.mask import Mask
from commons.layer import Layer
from commons.helpers import is_number
from commons.xml import *
from dataextractor import DataExtractor
from mapplot import mapplot

#Date extractor
def get_date_string(s):
    m = re.search('.*([0-9]{4})([0-9]{2})([0-9]{2}).*',s)
    longdate = "%s-%s-%s" % (m.group(1), m.group(2), m.group(3))
    shortdate = "%s%s%s" % (m.group(1), m.group(2), m.group(3))
    return longdate , shortdate

class Plot(object):
    def __init__(self, varname, layerlist, clim):
        #Input validation
        self.__varname = str(varname)
        if not isinstance(layerlist, (list, tuple)) or ((len(layerlist) > 0) and not isinstance(layerlist[0], (Layer,))):
            raise ValueError("layerlist must be a list of Layers")
        self.__layerlist = layerlist
        if not isinstance(clim, (list, tuple)) or (len(clim) != 2) or not (is_number(clim[0]) and is_number(clim[1])):
            raise ValueError("clim must be a list of two numbers")
        self.__clim = clim

    @property
    def varname(self):
        return self.__varname

    @property
    def layerlist(self):
        return self.__layerlist

    @property
    def clim(self):
        return self.__clim

    def append_layer(self, layer):
        if not isinstance(layer, (Layer,)):
            raise ValueError("layer must be a Layer object")
        self.__layerlist.append(layer)

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
                clim = literal_eval(get_node_attr(pdef, "clim"))
                plot = Plot(get_node_attr(pdef, "var"), [], clim)
                #For each depth element
                for d in get_subelements(pdef, "depth"):
                    plot.append_layer(Layer(get_node_attr(d, "top"), get_node_attr(d, "bottom")))
                self.__plotlist.append(plot)

    def plot_maps_data(self, min_ticks=4, max_ticks=8):
        for f in self.__netcdffileslist:
            longdate , shortdate = get_date_string(f)
            for p in self.__plotlist:
                de = DataExtractor(p.varname, f, self._mask)
                for l in p.layerlist:
                    outfilename = "%s/ave.%s.%s.%s.png" % (self.__outputdir,shortdate, p.varname, l)
                    mapdata = MapBuilder.get_layer_average(de, l)
                    fig = mapplot({'filename':f, 'varname':p.varname, 'clim':p.clim, 'layer':l, 'data':mapdata, 'date':longdate}, mask=self._mask, min_ticks=min_ticks, max_ticks=max_ticks)
                    fig.savefig(outfilename)
                    fig.clear()

    @staticmethod
    def get_layer_average(data_extractor, layer, timestep=0):
        """Returns a 2D NumPy array with the average weighted over depth.

        Args:
            - *data_extractor*: a DataExtractor object that contains the data to
              extract.
            - *layer*: a Layer object that contains the vertical boundaries for
              the computation.
            - *timestep* (optional): the index of the first dimension (time) in
              the model data. Default: 0.
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
            output = data_extractor.get_filled_values[timestep,top_index,:,:]
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
        v = np.array(data_extractor.values[timestep,top_index:bottom_index,:,:])
        #Build integral matrix (2D)
        integral = (v * dzm * lmask).sum(axis=0)
        #Build height matrix (2D)
        height = (dzm * lmask).sum(axis=0)
        indexmask = [height > 0]
        #Build output matrix (2D)
        output = np.full_like(integral, data_extractor.fill_value, dtype=np.double)
        output[indexmask] = integral[indexmask] / height[indexmask]
        return output

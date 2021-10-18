# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import os
import warnings
import numpy as np
from xml.dom import minidom
from ast import literal_eval

from commons.layer import Layer,LayerMap
from commons.utils import is_number, get_date_string
from commons.xml_module import *
from commons.dataextractor import DataExtractor
from commons.dataextractor import NotFoundError
from layer_integral.mapplot import mapplot,mapplot_medeaf_V5C,mapplot_nocolor
import datetime
import matplotlib.pyplot as pl

def warn_user(msg):
    warnings.warn(msg, SyntaxWarning, stacklevel=2)

class Plot(object):
    def __init__(self, varname, longvarname, units, layerlist, clim):
        #Input validation
        self.__varname = str(varname)
        self.__longvarname = str(longvarname)
        self.__units = units
        if not isinstance(layerlist, (list, tuple)) or ((len(layerlist) > 0) and not isinstance(layerlist[0], (Layer,))):
            raise ValueError("layerlist must be a list of Layers")
        self.__layerlist = layerlist
        if not (clim is None):
            if not isinstance(clim, (list, tuple)) or (len(clim) != 2) or not (is_number(clim[0]) and is_number(clim[1])):
                raise ValueError("clim must be a list of two numbers")
        self.__clim = clim
        self.__climslist = list()
        self.__depthfilters = list()

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
    @property
    def depthfilters(self):
        return self.__depthfilters

    def append_layer(self, layer, clim=None, mapdepthfilter=None):
        if not isinstance(layer, (Layer,)):
            raise ValueError("layer must be a Layer object")
        if mapdepthfilter is None:
            mapdepthfilter = layer.bottom
        self.__layerlist.append(layer)
        self.__climslist.append(clim)
        self.__depthfilters.append(mapdepthfilter)

class MapBuilder(object):

    def __init__(self, plotlistfile, TL, MaskObject, outputdir):
        '''
        Arguments : 
        * TL * is a Timelist object
        '''

        self.__TL =TL
        if os.path.isdir(outputdir):
            self.__outputdir = outputdir
        else:
            raise ValueError("outputdir must be a path to a directory")
        self._mask = MaskObject
        xmldoc = minidom.parse(plotlistfile)
        self.__plotlist = list()
        #For each LayersMaps element
        for lm in xmldoc.getElementsByTagName("LayersMaps"):
            #For each plots element
            for pdef in get_subelements(lm, "plots"):
                clim = get_node_attr(pdef, "clim")
                if not (clim is None):
                    clim = literal_eval(clim)
                plot = Plot(get_node_attr(pdef, "var"), get_node_attr(pdef, "longname"), get_node_attr(pdef, "plotunits"), [], clim)
                #For each depth element
                for d in get_subelements(pdef, "depth"):
                    clim = get_node_attr(d, "clim")
                    if not (clim is None):
                        clim = literal_eval(clim)
                    plot.append_layer(Layer(get_node_attr(d, "top"), get_node_attr(d, "bottom")), clim)
                self.__plotlist.append(plot)
    
    def read_background(self,filename):
        sfondo = pl.imread(filename)
        return sfondo
        
    def plot_maps_data(self, basemap_obj, maptype=0, background_img=None,nranks=1, rank=0):
        '''
        Parallel generator of a large set of images.
        Manages ave.date17.var.nc and ave.date17.phys.nc

        Arguments : 
          * basemap_obj *  a Basemap object
          * map_type      *  integer, flag used to choose a plot method in mapplot module
                = 0, the default, to call mapplot.mapplot()
                = 1, to call mapplot.mapplot_onlycolor()
                = 2, to call mapplot.mapplot_nocolor()
         * nranks, rank * integers to manage the ranks
        
        ''' 

        TL = self.__TL
        INPUTDIR = TL.inputdir
        nTimes   = TL.nTimes
        is_file_phys = TL.filelist[0].endswith("phys.nc")

        nplots = len(self.__plotlist)
        PROCESSES = np.arange(nTimes * nplots)
        for ip in PROCESSES[rank::nranks]:
            (itime, ivar) = divmod(ip,nplots)
            p = self.__plotlist[ivar]
            var     =  self.__plotlist[ivar].varname   #VARLIST[ivar]
            if is_file_phys:
                filename = INPUTDIR + "ave." + TL.Timelist[itime].strftime("%Y%m%d-%H:%M:%S.") + "phys.nc"
            else:
                filename= INPUTDIR + "ave." + TL.Timelist[itime].strftime("%Y%m%d-%H:%M:%S.") + var + ".nc"
            longdate , shortdate = get_date_string(filename)
#        for f in self.__netcdffileslist: 
#            for p in self.__plotlist:
            msg = "rank %d works on %s %s" %(rank,filename,var)
            print(msg)
            de = DataExtractor(self._mask, filename=filename, varname=p.varname)

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
                    fig, ax = mapplot_medeaf_V5C(mapdict, basemap_obj, self._mask, fig=None, ax=None , ncolors=24, logo=background_img)
                    fig.savefig(outfile + ".png",dpi=200)
                    pl.close(fig)

                if maptype == 2:
                    fig, ax = mapplot_nocolor({'varname':p.varname, 'clim':clim, 'layer':l, 'data':mapdata, 'date':longdate}, fig=fig, ax=ax, mask=self._mask, ncolors=24, coastline_lon=coastline_lon, coastline_lat=coastline_lat)
                    fig.canvas.print_figure(outfile + ".svg")
                    return


    @staticmethod
    def get_layer_max(data_extractor, layer):
        """Returns a 2D NumPy array with the maximum over depth.

        Args:
            - *data_extractor*: a DataExtractor object that contains the data to
              extract.
            - *layer*: a Layer object that contains the vertical boundaries for
              the computation.
        """
        if not isinstance(layer, (Layer,LayerMap)) :
            raise ValueError("layer must be a Layer or a MapLayer object")
        if not isinstance(data_extractor, (DataExtractor,)):
            raise ValueError("dataextractor must be a DataExtractor object")

        if isinstance(layer, (Layer,)) :
            #Find Z indices
            top_index = np.where(data_extractor._mask.zlevels >= layer.top)[0][0]
            bottom_index = np.where(data_extractor._mask.zlevels < layer.bottom)[0][-1]
            if top_index == bottom_index:
                #Just one layer so we return the sliced data
                output = data_extractor.filled_values[top_index,:,:]
                return output
            #Workaround for Python ranges
            bottom_index += 1
            #Get the slice of the values
            v = np.array(data_extractor.values[top_index:bottom_index,:,:])
            #Get the max of the values
            maximum = np.nanmax(v,0)
            output = maximum
            return output

        if isinstance(layer, (LayerMap,)) :
            #Find Z indices
            top_index    = np.zeros(layer.dimension,dtype=int)
            bottom_index = np.zeros(layer.dimension,dtype=int)
            for jj in range(layer.dimension[0]):
                for ii in range(layer.dimension[1]):
                    top_index[jj,ii] =  np.where(data_extractor._mask.zlevels >= layer.top[jj,ii])[0][0]
                    bottom_index[jj,ii] = np.where(data_extractor._mask.zlevels < layer.bottom[jj,ii])[0][-1]
                    #Workaround for Python ranges
                    bottom_index[jj,ii] += 1
            #Min top index
            min_top_index = top_index.min()
            #Max bottom index
            max_bottom_index = bottom_index.max()
            #Build local mask matrix
            lmask = np.array(data_extractor._mask.mask[min_top_index:max_bottom_index,:,:], dtype=np.double)
            for jj in range(layer.dimension[0]):
                for ii in range(layer.dimension[1]):
                    j = top_index[jj,ii] - min_top_index 
                    lmask[:j,jj,ii] = 0
                    i = bottom_index[jj,ii] - min_top_index
                    lmask[i:,jj,ii] = 0 
            #Get the slice of the values
            v = np.array(data_extractor.values[min_top_index:max_bottom_index,:,:])
            v[lmask==0] = 0
            #Build max matrix (2D)
            maximum = np.nanmax(v,0)
            output = maximum
            return output

    @staticmethod
    def get_layer_average(data_extractor, layer):
        """Returns a 2D NumPy array with the average weighted over depth.

        Args:
            - *data_extractor*: a DataExtractor object that contains the data to
              extract.
            - *layer*: a Layer object that contains the vertical boundaries for
              the computation.
        """
        if not isinstance(layer, (Layer,LayerMap)):
            raise ValueError("layer must be a Layer object")
        if not isinstance(data_extractor, (DataExtractor,)):
            raise ValueError("dataextractor must be a DataExtractor object")
        if isinstance(layer, (Layer,)) :
            #Find Z indices
            top_index = np.where(data_extractor._mask.zlevels >= layer.top)[0][0]
            if layer.bottom < data_extractor._mask.zlevels[0]:
                bottom_index=0
            else:
                bottom_index = np.where(data_extractor._mask.zlevels < layer.bottom)[0][-1]
            if bottom_index == (top_index-1) : #when layer.top==layer.bottom
                bottom_index = top_index
            if top_index == bottom_index:
                #Just one layer so we return the sliced data
                output = data_extractor.filled_values[top_index,:,:]
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
            indexmask = height > 0
            #Build output matrix (2D)
            output = np.full_like(integral, data_extractor.fill_value, dtype=np.double)
            output[indexmask] = integral[indexmask] / height[indexmask]
            return output

        if isinstance(layer, (LayerMap,)):
            #Find Z indices
            top_index    = np.zeros(layer.dimension,dtype=int)
            bottom_index = np.zeros(layer.dimension,dtype=int)
            for jj in range(layer.dimension[0]):
                for ii in range(layer.dimension[1]):
                    top_index[jj,ii] =  np.where(data_extractor._mask.zlevels >= layer.top[jj,ii])[0][0]
                    bottom_index[jj,ii] = np.where(data_extractor._mask.zlevels < layer.bottom[jj,ii])[0][-1]
                    #Workaround for Python ranges
                    bottom_index[jj,ii] += 1
            #Min top index
            min_top_index = top_index.min()
            #Max bottom index
            max_bottom_index = bottom_index.max()
            #Build local mask matrix
            lmask = np.array(data_extractor._mask.mask[min_top_index:max_bottom_index,:,:], dtype=np.double)
            for jj in range(layer.dimension[0]):
                for ii in range(layer.dimension[1]):
                    j = top_index[jj,ii] - min_top_index 
                    lmask[:j,jj,ii] = 0
                    i = bottom_index[jj,ii] - min_top_index
                    lmask[i:,jj,ii] = 0 
            #Build dz matrix
            dzm = np.ones_like(lmask, dtype=np.double)
            j = 0
            for i in range(min_top_index, max_bottom_index):
                dzm[j,:,:] = data_extractor._mask.dz[i]
                j += 1
            #Get the slice of the values
            v = np.array(data_extractor.values[min_top_index:max_bottom_index,:,:])
            #Build integral matrix (2D)
            integral = (v * dzm * lmask).sum(axis=0)
            #Build height matrix (2D)
            height = (dzm * lmask).sum(axis=0)
            indexmask = height > 0
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
        if not isinstance(layer, (Layer,LayerMap)):
            raise ValueError("layer must be a Layer object")
        if not isinstance(data_extractor, (DataExtractor,)):
            raise ValueError("dataextractor must be a DataExtractor object")
        #Find Z indices
        if isinstance(layer, (Layer,)):
            top_index = np.where(data_extractor._mask.zlevels >= layer.top)[0][0]
            if layer.bottom < data_extractor._mask.zlevels[0]:
                bottom_index=0
            else:
                bottom_index = np.where(data_extractor._mask.zlevels < layer.bottom)[0][-1]
            if top_index == bottom_index:
                #Just one layer so we return the sliced data
                output = data_extractor.filled_values[top_index,:,:]
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
            indexmask = height > 0
            #Build output matrix (2D)
            output = np.full_like(integral, data_extractor.fill_value, dtype=np.double)
            output[indexmask] = integral[indexmask]
            return output

        if isinstance(layer, (LayerMap,)):
           #Find Z indices
            top_index    = np.zeros(layer.dimension,dtype=int)
            bottom_index = np.zeros(layer.dimension,dtype=int)
            for jj in range(layer.dimension[0]):
                for ii in range(layer.dimension[1]):
                    top_index[jj,ii] =  np.where(data_extractor._mask.zlevels >= layer.top[jj,ii])[0][0]
                    bottom_index[jj,ii] = np.where(data_extractor._mask.zlevels < layer.bottom[jj,ii])[0][-1]
                    #Workaround for Python ranges
                    bottom_index[jj,ii] += 1
            #Min top index
            min_top_index = top_index.min()
            #Max bottom index
            max_bottom_index = bottom_index.max()
            #Build local mask matrix
            lmask = np.array(data_extractor._mask.mask[min_top_index:max_bottom_index,:,:], dtype=np.double)
            for jj in range(layer.dimension[0]):
                for ii in range(layer.dimension[1]):
                    j = top_index[jj,ii] - min_top_index
                    lmask[:j,jj,ii] = 0
                    i = bottom_index[jj,ii] - min_top_index
                    lmask[i:,jj,ii] = 0
            #Build dz matrix
            dzm = np.ones_like(lmask, dtype=np.double)
            j = 0
            for i in range(min_top_index, max_bottom_index):
                dzm[j,:,:] = data_extractor._mask.dz[i]
                j += 1
            #Get the slice of the values
            v = np.array(data_extractor.values[min_top_index:max_bottom_index,:,:])
            #Build integral matrix (2D)
            integral = (v * dzm * lmask).sum(axis=0)
            #Build height matrix (2D)
            height = (dzm * lmask).sum(axis=0)
            indexmask = [height > 0]
            #Build output matrix (2D)
            output = np.full_like(integral, data_extractor.fill_value, dtype=np.double)
            output[indexmask] = integral[indexmask]
            return output

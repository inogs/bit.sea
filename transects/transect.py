# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import warnings
import numpy as np
from xml.dom import minidom
from ast import literal_eval
from commons.segment import Segment
from commons.utils import is_number
from commons.xml_module import *
from commons.mask import Mask
from commons.dataextractor import DataExtractor

class Transect(object):
    """Stores a multiple segment transect definition.
    """
    def __init__(self, varname, clim, segmentlist, name=None, unit=None):
        """Transect constructor.

        Args:
            - *varname*: the name of the variable to extract.
            - *clim*: a list or tuple with the two limit values (minimum and
              maximum).
            - *segmentlist*: a list Segment objects. Can be an empty list.
            - *name* (optional): a string defining the transect's name. If set
              to None it will be set to varname (default: None).
            - *unit* (optional): a string defining the transect's measurement
              unit (default: None).
        """
        if varname is None:
            raise ValueError("varname cannot be None")
        self.__varname = str(varname)
        if name is None:
            self.__name = self.__varname
        else:
            self.__name = str(name)
        if unit is None:
            self.__unit = None
        else:
            self.__unit = str(unit)
        if not isinstance(clim, (list, tuple)) or (len(clim) != 2) or not (is_number(clim[0]) and is_number(clim[1])):
            raise ValueError("clim must be a list of two numbers")
        self.__clim = clim
        if not isinstance(segmentlist, (list, tuple)) or ((len(segmentlist) > 0) and not isinstance(segmentlist[0], (Segment,))):
            raise ValueError("segmentlist must be a list of Segments")
        self.__segmentlist = segmentlist
        #Variable data cache
        self.__datacache = { 'filename':None, 'data':None }
        self.__mask = None

    @property
    def name(self):
        return self.__name

    @property
    def varname(self):
        return self.__varname

    @property
    def clim(self):
        return self.__clim

    @property
    def unit(self):
        return self.__unit

    @property
    def segmentlist(self):
        return self.__segmentlist

    def fill_data_cache(self, mask, data):
        """Fills the object cache with data from a Numpy array

        Args:
            - *mask*: a commons.Mask object.
            - *data*: a numpy.ndarray containing the data to slice.
        """
        if not isinstance(mask, (Mask,)):
            raise ValueError("mask must be a Mask object.")
        self.__mask = mask
        if not isinstance(data, np.ndarray):
            raise ValueError("data must be a Numpy array")
        if len(data.shape) != 3:
            raise ValueError("data must be a tridimensional matrix")
        self.__datacache['filename'] = None
        self.__datacache['data'] = data

    def fill_data_cache_from_file(self, mask, datafilepath, fill_value=np.nan):
        """Fills the object cache with data from a NetCDF file

        Args:
            - *mask*: a commons.Mask object.
            - *datafilepath*: path to the NetCDF data file.
            - *fill_value* (optional): value to use when there's no data
              available (default: np.nan).
        """
        if not isinstance(mask, (Mask,)):
            raise ValueError("mask must be a Mask object.")
        self.__mask = mask
        #Check if we don't have cached data
        datafilepath = str(datafilepath)
        if (self.__datacache['filename'] is None) or (self.__datacache['filename'] != datafilepath):
            #Read data from the NetCDF file
            de = DataExtractor(mask, filename=datafilepath, varname=self.__varname, fill_value=fill_value)
            self.__datacache['filename'] = datafilepath
            self.__datacache['data'] = de.filled_values

    def get_segment_data(self):
        """Retrieves the requested segments data.

        Returns: a list of dictionaries. Each dictionary contains the following keys:
                 - segment: the Segment instance
                 - data: the corresponding data (Numpy ndarray).
                 - h_vals: the degrees of each point.
                 - z_vals: the depth in meters of each point.
        """
        if self.__mask is None:
            raise ValueError("Missing mask data. Try call fill_data_cache before asking for segment data")
        output = list()
        for i,s in enumerate(self.__segmentlist):
            data, h_vals, z_vals = self.__get_segment_data(i)
            output.append({'segment': s, 'data': data, 'h_vals': h_vals, 'z_vals': z_vals})
        return output

    def get_transect_data(self, segment_index, mask, datafilepath, fill_value=np.nan):
        """Extracts transect data from a NetCDF file.

        Args:
            - *segment_index*: the index of the segment to retrieve.
            - *mask*: a commons.Mask object.
            - *datafilepath*: path to the NetCDF data file.
            - *fill_value* (optional): value to use when there's no data
              available (default: np.nan).
            - *timestep* (optional): the time index (default: 0).
        Returns: a NumPy array that contains the requested data.
        """
        self.fill_data_cache_from_file(mask, datafilepath, fill_value)
        data, _, _ = self.__get_segment_data(segment_index)
        return data

    @staticmethod
    def get_list_from_XML_file(plotlistfile):
        xmldoc = minidom.parse(plotlistfile)
        output = list()
        #For each Transects element
        for t in xmldoc.getElementsByTagName("Transects"):
            #Build the segment list
            segmentlist = list()
            for sdef in get_subelements(t, "transect"):
                sname = get_node_attr(sdef, "name")
                lon_min = literal_eval(get_node_attr(sdef, "lonmin"))
                lat_min = literal_eval(get_node_attr(sdef, "latmin"))
                lon_max = literal_eval(get_node_attr(sdef, "lonmax"))
                lat_max = literal_eval(get_node_attr(sdef, "latmax"))
                segmentlist.append(Segment((lon_min, lat_min),(lon_max, lat_max),sname))
            #For each vars list
            for vl in get_subelements(t, "vars"):
                for v in get_subelements(vl, "var"):
                    varname = get_node_attr(v, "name")
                    if varname is None:
                        raise SyntaxError("name is not defined")
                    clim = get_node_attr(v, "clim")
                    if clim is None:
                        raise ValueError("clim not set for variable %s" % (varname,))
                    try:
                        clim = literal_eval(get_node_attr(v, "clim"))
                    except SyntaxError as e:
                        msg = "Variable %s invalid syntax for clim: %s" % (varname, e.text)
                        raise SyntaxError(msg)
                    except:
                        raise
                    output.append(Transect(varname, clim, segmentlist))
        return output

    #Private methods
    def __get_segment_data(self, segment_index):
        #Input validation
        if not is_number(segment_index):
            raise ValueError("'%s' is not a number." % (segment_index,))
        seg = self.__segmentlist[segment_index]
        #Retrieve indices
        x_min, y_min = self.__mask.convert_lon_lat_to_indices(seg.lon_min, seg.lat_min)
        x_max, y_max = self.__mask.convert_lon_lat_to_indices(seg.lon_max, seg.lat_max)
        #Check for single point case
        if (x_min == x_max) and (y_min == y_max):
            raise ValueError("Invalid segment: %s" % (seg,))
        #Get the output data
        if x_min == x_max:
            data = np.array(self.__datacache['data'][:, y_min:y_max, x_min], dtype=float)
            h_vals = self.__mask.ylevels[y_min:y_max,x_min]
        elif y_min == y_max:
            data = np.array(self.__datacache['data'][:, y_min, x_min:x_max], dtype=float)
            h_vals = self.__mask.xlevels[y_min,x_min:x_max]
        else:
            #Divide the segment in points
            xs = np.linspace(x_min,x_max,num=seg.points)
            ys = np.linspace(y_min,y_max,num=seg.points)
            xs = np.round(xs).astype(np.int32)
            ys = np.round(ys).astype(np.int32)
            #Fill data with the nearest neighbour
            data = list()
            #print "seg.points = ", seg.points
            for i in range(seg.points):
                #print "indexes=", i, ys[i], xs[i]
                data.append(self.__datacache['data'][:,ys[i],xs[i]])
            #Convert data to a Numpy array
            data = np.array(data)
            data = data.transpose()
            #TODO: find a meaningful array for h_vals.
            h_vals = np.array(range(seg.points))
        z_vals = self.__mask.zlevels
        return data, h_vals, z_vals


# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
from xml.dom import minidom
from ast import literal_eval
from commons.segment import Segment
from commons.helpers import is_number
from commons.xml import *

class Transect(object):
    """Stores a multiple segment transect definition.
    """
    def __init__(self, varname, clim, segmentlist):
        """Transect constructor.

        Args:
            - *varname*: the name of the variable to extract.
            - *clim*: a list or tuple with the two limit values (minimum and
              maximum)
            - *segmentlist*: a list Segment objects. Can be an empty list.
        """
        self.__varname = str(varname)
        if not isinstance(clim, (list, tuple)) or (len(clim) != 2) or not (is_number(clim[0]) and is_number(clim[1])):
            raise ValueError("clim must be a list of two numbers")
        self.__clim = clim
        if not isinstance(segmentlist, (list, tuple)) or ((len(segmentlist) > 0) and not isinstance(segmentlist[0], (Segment,))):
            raise ValueError("segmentlist must be a list of Segments")
        self.__segmentlist = segmentlist

    @property
    def varname(self):
        return self.__varname

    @property
    def clim(self):
        return self.__clim

    @property
    def segmentlist(self):
        return self.__segmentlist

    @staticmethod
    def get_list_from_XML_file(plotlistfile):
        xmldoc = minidom.parse(plotlistfile)
        output = list()
        #For each Transects element
        for t in xmldoc.getElementsByTagName("Transects"):
            #Build the segment list
            segmentlist = list()
            for sdef in get_subelements(t, "transect"):
                lon_min = literal_eval(get_node_attr(sdef, "lonmin"))
                lat_min = literal_eval(get_node_attr(sdef, "latmin"))
                lon_max = literal_eval(get_node_attr(sdef, "lonmax"))
                lat_max = literal_eval(get_node_attr(sdef, "latmax"))
                if (lon_min != lon_max) and (lat_min != lat_max):
                    raise NotImplementedError("You have to fix a coordinate. Either Longitude or Latitude.")
                segmentlist.append(Segment((lon_min, lat_min),(lon_max, lat_max)))
            #For each vars list
            for vl in get_subelements(t, "vars"):
                for v in get_subelements(vl, "var"):
                    varname = get_node_attr(v, "name")
                    clim = literal_eval(get_node_attr(v, "clim"))
                    output.append(Transect(varname, clim, segmentlist))
        return output

# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np
class Layer(object):
    def __init__(self,top, bottom):
        t = float(top)
        b = float(bottom)
        if t <= b:
            self.__top = t
            self.__bottom = b
        else:
            raise ValueError("top must be above of bottom")

    def __repr__(self):
        if self.top ==self.bottom:
            return "Layer %g m" %self.__top
        return "Layer %g-%g m" %(self.__top, self.__bottom)

    def __str__(self):
        if self.top ==self.bottom:
            return "Layer %g m" %self.__top
        return "layer%g-%g" %(self.__top, self.__bottom)

    def string(self):
        if self.top ==self.bottom:
            return "%gm" %self.__top
        return "%g-%gm" %(self.__top, self.__bottom)

    def longname(self):
        if self.top ==self.bottom:
            return "%04gm" %self.__top
        return "%04g-%04gm" %(self.__top, self.__bottom)

    @property
    def top(self):
        return self.__top

    @property
    def bottom(self):
        return self.__bottom

class LayerMap(object):
    def __init__(self, mask, top, bottom):
        if (top.shape==mask.shape[1:]) & (bottom.shape==mask.shape[1:]):
            self.__mask = mask
            self.__dim = mask.shape[1:]
        else:
            raise ValueError("top and bottom dimensions must be equal to mask dimensions along lat and lon")

        if np.all(top <= bottom):
            self.__top = top
            self.__bottom = bottom
        else:
            raise ValueError("top must be above of bottom")

    def __repr__(self):
        return "Map of Layers with dimensions %g,%g" %(self.__dim[0], self__dim[1])

    def __str__(self):
        return "maplayer(%g,%g)" %(self.__dim[0], self__dim[1])

    def string(self):
        return "maplayer(%g,%g)" %(self.__dim[0], self__dim[1])

    def longname(self):
        return "Map of Layers, dimension %g,%g" %(self.__dim[0], self__dim[1])

    @property
    def top(self):
        return self.__top

    @property
    def bottom(self):
        return self.__bottom

    @property
    def mask(self):
        return self.__mask

    @property
    def dimension(self):
        return self.__dim


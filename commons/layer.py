# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
class Layer(object):
    def __init__(self,top, bottom):
        t = float(top)
        b = float(bottom)
        if t < b:
            self.__top = t
            self.__bottom = b
        else:
            raise ValueError("top must be above of bottom")

    def __repr__(self):
        return "Layer %g-%g m" %(self.__top, self.__bottom)

    def __str__(self):
        return "layer%g-%g" %(self.__top, self.__bottom)

    @property
    def top(self):
        return self.__top

    @property
    def bottom(self):
        return self.__bottom


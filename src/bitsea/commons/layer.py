# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
from numbers import Real

import numpy as np


class Layer:
    """
    This object represents a depth interval in meters. It is defined
    by two numbers, `top` and `bottom`
    """

    def __init__(self, top: Real, bottom: Real):
        t = float(top)
        b = float(bottom)
        if t > b:
            raise ValueError("top must be above of bottom")

        self.__top = t
        self.__bottom = b

    def __repr__(self):
        return f"Layer({self.__top}, {self.__bottom})"

    def __str__(self):
        if self.top == self.bottom:
            return f"Layer {self.__top} m"
        return f"Layer {self.__top:g}-{self.__bottom:g} m"

    def __eq__(self, other):
        if not isinstance(other, Layer):
            return NotImplemented
        top_close = np.abs(self.top - other.top) < 1e-6
        bottom_close = np.abs(self.bottom - other.bottom) < 1e-6
        return top_close and bottom_close

    def string(self):
        if self.top == self.bottom:
            return f"{self.__top:g}m"
        return f"{self.__top:g}-{self.__bottom:g}m"

    def longname(self):
        if self.top == self.bottom:
            return f"{self.__top:04g}m"
        return f"{self.__top:04g}-{self.__bottom:04g}m"

    def to_file_postfix(self):
        """
        Returns a filename-safe string representation of this layer.

        Returns a string that describes the layer in a format suitable for use
        in filenames, containing only alphanumeric characters, dots, and
        hyphens.
        """
        if self.top == self.bottom:
            return f"layer{self.__top:g}"
        return f"layer{self.__top:g}-{self.__bottom:g}"

    @property
    def top(self):
        return self.__top

    @property
    def bottom(self):
        return self.__bottom


class LayerMap(object):
    def __init__(self, mask, top, bottom):
        if (top.shape == mask.shape[1:]) & (bottom.shape == mask.shape[1:]):
            self.__mask = mask
            self.__dim = mask.shape[1:]
        else:
            raise ValueError(
                "top and bottom dimensions must be equal to mask dimensions along lat and lon"
            )

        if np.all(top <= bottom):
            self.__top = top
            self.__bottom = bottom
        else:
            raise ValueError("top must be above of bottom")

    def __repr__(self):
        return "Map of Layers with dimensions %g,%g" % (
            self.__dim[0],
            self.__dim[1],
        )

    def __str__(self):
        return "maplayer(%g,%g)" % (self.__dim[0], self.__dim[1])

    def string(self):
        return "maplayer(%g,%g)" % (self.__dim[0], self.__dim[1])

    def longname(self):
        return "Map of Layers, dimension %g,%g" % (self.__dim[0], self.__dim[1])

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

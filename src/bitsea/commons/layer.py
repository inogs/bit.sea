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
        return f"{self.__class__.__name__}({self.__top}, {self.__bottom})"

    def __str__(self):
        if self.top == self.bottom:
            return f"{self.__class__.__name__} {self.__top} m"
        return f"{self.__class__.__name__} {self.__top:g}-{self.__bottom:g} m"

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

    @property
    def top(self):
        return self.__top

    @property
    def bottom(self):
        return self.__bottom


class LayerMap(Layer):
    def __init__(self, mask, top: Real, bottom: Real):
        super().__init__(top, bottom)
        self.__mask = mask

    @property
    def mask(self):
        return self.__mask

    @property
    def dimension(self):
        return self.__mask.shape[1:]

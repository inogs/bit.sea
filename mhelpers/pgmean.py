# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np
from scipy.signal import gaussian
from mean import GaussianMean
from math import exp

class PGaussianMean(GaussianMean):
    '''
    Gaussian weighted against pressure moving average helper object
    '''

    def compute(self, values, pressure_values):
        if self._check_compute_input(values, pressure_values, pressure_required=True) == False:
            return np.array(values)
        l = len(values)
        #Ensure we have np.arrays
        values = np.array(values, dtype=float)
        pressure_values = np.array(pressure_values, dtype=float)
        output = np.empty((l,), dtype=float)
        #Build base gaussian weights
        w = gaussian(self._i, self._sigma)
        for i in range(l):
            rbegin = (i - (self._i // 2))
            rend = (i + (self._i // 2)) + 1
            wbegin = 0
            if rbegin < 0:
                wbegin = -rbegin
                rbegin = 0
            if rend > l:
                rend = l
            wend = (rend - rbegin) + wbegin
            #Build pressure-modified weights
            wm = w.copy()
            pvals = pressure_values[rbegin:rend]
            k = 0
            #Build distances vector
            dvals = np.absolute(pressure_values[rbegin:rend] - pressure_values[i])
            dmin = np.inf
            for d in dvals:
                if d > 0 and d < dmin:
                   dmin = d
            for j in range(wbegin, wend):
                d = dvals[k]
                if d == 0:
                    d = 1
                else:
                    d = (dmin / d)
                wm[j] = wm[j] * d
                k = k + 1
            output[i] = sum((values[rbegin:rend] * wm[wbegin:wend])) / sum(wm[wbegin:wend])
        return output

class PLGaussianMean(GaussianMean):

    def _gaussian(self, value):
        return exp(-0.5*(value / self._sigma)**2)

    def compute(self, values, pressure_values):
        if self._check_compute_input(values, pressure_values, pressure_required=True) == False:
            return np.array(values)
        l = len(values)
        #Ensure we have np.arrays
        values = np.array(values, dtype=float)
        pressure_values = np.array(pressure_values, dtype=float)
        output = np.empty((l,), dtype=float)
        w = np.zeros(self._i, dtype=float)
        for i in range(l):
            #Build indices
            rbegin = (i - (self._i // 2))
            rend = (i + (self._i // 2)) + 1
            wbegin = 0
            if rbegin < 0:
                wbegin = -rbegin
                rbegin = 0
            if rend > l:
                rend = l
            wend = (rend - rbegin) + wbegin
            #Build distances vector
            dvals = np.absolute(pressure_values[rbegin:rend] - pressure_values[i])
            dmax = max(dvals)
            k = 0
            #Build pressure-modified weights
            for j in range(wbegin, wend):
                w[j] = self._gaussian(dvals[k]/dmax)
                k = k + 1
            #Perform the weighted average
            output[i] = sum((values[rbegin:rend] * w[wbegin:wend])) / sum(w[wbegin:wend])
        return output

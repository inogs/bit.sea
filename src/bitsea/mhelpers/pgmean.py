# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import numpy as np
from scipy.signal import gaussian
from mhelpers.mean import GaussianMean
from math import exp

class PGaussianMean(GaussianMean):
    """Gaussian weighted against pressure moving average helper object

    This class extends GaussianMean and modifies the way the weights are
    computed taking into account the pressure values. First the weights are
    computed like the ones used in GaussianMean, then a distance vector is
    computed for the samples in the current interval. The distances are then
    normalized according to the following formula::

        if d[i] != 0:
            d[i] = d_min / d[i]
        else:
            d[i] = 1

    Where d_min is the minimum non-zero distance. E.G.::

        [3, 2, 1, 0, 10, 20, 30]

    becomes::

        [1/3, 1/2, 1, 1, 1/10, 1/20, 1/30]

    This new vector is then multiplied element-by-element with the weight
    vector and these new weights are then used to compute the moving average::

            output[i] = sum((values[(i - interval // 2):(i + (interval // 2 ) + 1)] * weights)) / sum(weights)

    """

    def compute(self, values, pressure_values):
        """Performs the computation.

        Args:
            - *values*: the array of values to smooth.
            - *pressure_values*: the array of pressure values to be used on the
              weights. MUST have the same length of values.
        """
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
    """Logarithmic-like smoothing helper object.

    This class smooths the values according to the following algorithm:

    1. Build a distance vector as long as the current interval based on
       pressure_values and centered on the current sample.
    2. Normalize the distances::

        d[i] = d[i] / d_max

       where d_max is the maximum distance.
    3. Compute the weights based on the normalized distances::

        weights[i] = exp(-0.5*(d[i] / sigma)**2)

    4. Compute the output::

            output[i] = sum((values[(i - interval // 2):(i + (interval // 2 ) + 1)] * weights)) / sum(weights)

    Basically instead of weighting the weights it takes different points on the
    Gauss bell curve according to the distance from the current sample.
    """

    def _gaussian(self, value):
        return exp(-0.5*(value / self._sigma)**2)

    def compute(self, values, pressure_values):
        """Performs the computation.

        Args:
            - *values*: the array of values to smooth.
            - *pressure_values*: the array of pressure values to be used on the
              weights. MUST have the same length of values.
        """
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
            #Prevent division by zero error
            if dmax == 0:
                dmax = 1.0
            k = 0
            #Build pressure-modified weights
            for j in range(wbegin, wend):
                w[j] = self._gaussian(dvals[k]/dmax)
                k = k + 1
            #Perform the weighted average
            output[i] = sum((values[rbegin:rend] * w[wbegin:wend])) / sum(w[wbegin:wend])
        return output

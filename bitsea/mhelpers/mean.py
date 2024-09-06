# Copyright (c) 2015 eXact Lab srl
# Authors: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
#          Stefano Piani <stefano.piani@exact-lab.it>
import numpy as np
from scipy.signal import gaussian

class Mean(object):
    """Base class for data filtering.

    This is an abstract class, do not attempt to create a Mean object.
    """
    def __init__(self, interval):
        """Mean class constructor.

        Args:
            - *interval*: the amount of nearby points that the filter has to
              consider for each point of the values' array including the point
              itself.

        Raises:
            - *NotImplementedError* because this is an abstract class.
        """
        raise NotImplementedError()

    def compute(self, values, pressure_values):
        """Performs the filtering.

        Args:
            - *values*: NumPy array, list or tuple of values. Must have the
              same length as pressure_values.
            - *pressure_values*: NumPy array, list or tuple of pressure values.
              Must have the same length as values.

        Raises:
            - *NotImplementedError* because this is an abstract class.
        """
        raise NotImplementedError()

    def _check_compute_input(self, values, pressure_values, pressure_required=False):
        """Helper function to perform input validation and avoid meaningless computation.

        Args:
            - *values*: the same input vector that has been passed to compute.
            - *pressure_values*: the same input vector that has been passed to compute.
            - *pressure_required*: if False pressure_values can be set to None.

        Raises:
            - *TypeError* if values is not an array of numbers.
        """
        assert not values is None
        if pressure_required:
            assert (not (pressure_values is None)) and (len(values) == len(pressure_values))
        else:
            assert (pressure_values is None) or (len(values) == len(pressure_values))
        l = len(values)
        if l==0:
            return False
        if not isinstance(values[0], (int, float, complex, np.float32)):
            raise TypeError()
        if l == 1 or self._i == 0:
            return False
        return True

class GaussianMean(Mean):
    """Gaussian weighted moving average helper object.
    """
    def __init__(self, interval, sigma=1.0):
        """GaussianMean class constructor.

        Args:
            - *interval*: the amount of nearby points that the filter has to
              consider for each point of the values' array including the point
              itself. The interval is centered on the current point and an even
              interval will be promoted to the next odd value (e.g. 4 will
              become 5 so two points before and two points after plus the
              current point will be considered).
            - *sigma*: the sigma parameter for the Gauss bell curve.
              Defaults to 1.0 .

        Raises:
            - *ValueError* if interval or sigma are negative or if they cannot
              be converted to a number.
        """
        if isinstance(interval, int ):
            self._i = interval
        elif isinstance(int(interval), int):
            self._i = int(interval)
        if self._i < 0:
            raise ValueError("interval should be positive")
        if isinstance(float(sigma), (float,)):
            sigma = float(sigma)
            if sigma >=0:
                self._sigma = sigma
            else:
                raise ValueError("sigma should be positive")

    def compute(self, values, pressure_values=None):
        """Performs the gaussian filtering based on values' indices.

        The filtering algorithm uses a moving average weighted according to this
        formula::

            weights[i] = exp(-0.5*((i - interval // 2) / sigma)**2)

        For a three-point interval with sigma of 1.0 this leads to the following
        weight vector::

            [ 0.60653066,  1.        ,  0.60653066]

        The output vector is then computed as::

            output[i] = sum((values[(i - interval // 2):(i + (interval // 2 ) + 1)] * weights)) / sum(weights)

        Args:
            - *values*: the input array of values to filter.
            - *pressure_values*: ignored (can be set to None). If it is set it
              MUST have the same lenght of values' array.
        """
        if self._check_compute_input(values, pressure_values) == False:
            return np.array(values)
        l = len(values)
        #Ensure we have np.arrays
        values = np.array(values, dtype=float)
        output = np.empty((l,), dtype=float)
        #Build gaussian weights
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
            output[i] = sum((values[rbegin:rend] * w[wbegin:wend])) / sum(w[wbegin:wend])
        return output

class MovingAverage(Mean):
    """Simple in-place moving average.

    This class implements the classic moving average algorithm.
    Each sample will be summed with all the samples in the interval and the
    sum will be divided by the length of the interval, then the result will
    replace the value of the sample and the algorithm will evaluate the next
    interval.
    """

    def __init__(self, interval):
        """MovingAverage class constructor.

        Args:
            - *interval*: the interval of samples to consider.

        Raises:
            - ValueError if interval is negative.
        """
        if isinstance(interval, int):
            self._i = interval
        elif isinstance(int(interval), int):
            self._i = int(interval)
        if self._i < 0:
            raise ValueError("interval should be positive")

    def compute(self, values, pressure_values=None):
        """Returns the moving average on values.

        Args:
            - *values*: the input array of values to filter.
            - *pressure_values*: ignored (can be set to None). If it is set it
              MUST have the same lenght of values' array.
        """
        if self._check_compute_input(values, pressure_values) == False:
            return np.array(values)
        l = len(values)
        #Ensure we have np.arrays
        values = np.array(values, dtype=float)
        output = values.copy()
        for i in range(l):
            rbegin = (i - (self._i // 2))
            rend = (i + (self._i // 2)) + 1
            if rbegin < 0:
                rbegin = 0
            if rend > l:
                rend = l
            output[i] = sum(output[rbegin:rend]) / (rend-rbegin)
        return output

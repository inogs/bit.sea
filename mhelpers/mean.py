import numpy as np
from scipy.signal import gaussian

class Mean(object):
    def __init__(self, interval):
        raise NotImplementedError()

    def compute(self, values, pressure_values):
        raise NotImplementedError()

class GaussianMean(Mean):
    '''
    Gaussian weighted moving average helper object
    '''
    def __init__(self, interval, sigma=1.0):
        if isinstance(interval, (int, long )):
            self._i = interval
        elif isinstance(int(interval), (int, long)):
            self._i = int(interval)
        if self._i < 0:
            raise ValueError("interval should be positive")
        if isinstance(float(sigma), (float,)):
            sigma = float(sigma)
            if sigma >=0:
                self._sigma = sigma
            else:
                raise ValueError("sigma should be positive")

    def compute(self, values, pressure_values):
        l = len(values)
        if l==0:
            return np.array(values)
        if not isinstance(values[0], (int, long, float, complex)):
            raise TypeError()
        if l == 1 or self._i == 0:
            return np.array(values)
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

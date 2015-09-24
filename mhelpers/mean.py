import numpy as np
from scipy.signal import gaussian

'''
Moving average helper object
'''
class mean():
    def __init__(self, interval, sigma=1.0):
        if isinstance(interval, (int, long )):
            self._i = interval
        elif isinstance(int(interval), (int, long)):
            self._i = int(interval)
        if self._i < 0:
            raise ValueError("interval should be positive")
        if isinstance(sigma, (int, long, float)) or isinstance(float(sigma), (float,)):
            if sigma >=0:
                self._sigma = sigma
            else:
                raise ValueError("sigma should be positive")
        else:
            raise ValueError("sigma should be a float")

    def compute(self, steps, values):
        if not len(steps) == len(values):
            raise ValueError("steps and values arrays must have the same length")
        l = len(steps)
        if not isinstance(steps[0], (int, long, float, complex)):
            raise TypeError()
        if not isinstance(values[0], (int, long, float, complex)):
            raise TypeError()
        if l == 1 or self._i == 0:
            return values
        #Ensure we have np.arrays
        steps = np.array(steps, dtype=float)
        values = np.array(values, dtype=float)
        output = np.empty([l,], dtype=float)
        for i in range(l):
            rbegin = (i - (self._i/2))
            rend = (i + (self._i/2)) + 1
            if rbegin < 0:
                rbegin = 0
            if rend > l:
                rend = l
            #Build gaussian weights
            w = gaussian((rend-rbegin), self._sigma)
            output[i] = sum((values[rbegin:rend] * w)) / sum(w)
        return output

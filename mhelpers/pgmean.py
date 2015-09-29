import numpy as np
from scipy.signal import gaussian
from mean import GaussianMean

class PGaussianMean(GaussianMean):
    '''
    Gaussian weighted against pressure moving average helper object
    '''

    def compute(self, values, pressure_values):
        if self._check_compute_input(values, pressure_values) == False:
            return np.array(values)
        if pressure_values == None:
            raise TypeError("pressure_values must be defined")
        l = len(values)
        #Ensure we have np.arrays
        values = np.array(values, dtype=float)
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
            for j in range(wbegin, wend):
                d = abs(pvals[k] - pressure_values[i])
                if d == 0:
                    d = 1
                wm[j] = wm[j] / d
                k = k + 1
            print "wm: " + str(wm[wbegin:wend])
            output[i] = sum((values[rbegin:rend] * wm[wbegin:wend])) / sum(wm[wbegin:wend])
        return output

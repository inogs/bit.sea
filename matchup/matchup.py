import numpy as np
import pylab as pl

class Layer(object):
    def __init__(self, top, bottom):
        self.__top = top
        self.__bottom = bottom
    def __str__(self):
        return "Layer %d-%d m" %(self.__top, self.__bottom)

    @property
    def top(self):
        return self.__top

    @property
    def bottom(self):
        return self.__bottom

class matchup():
    def __init__(self, Model, Ref):
        self.Model = Model
        self.Ref   = Ref

    def diff(self):
        return self.Ref - self.Model

    def number(self):        
        return len(self.Model)

    def variances(self):
        return self.Ref.std()**2, self.Model.std()**2

    def medians(self):
        '''
        return median(Reference), median(Model)
        '''
        return np.median(self.Ref), np.median(self.Model)

    def covariance(self):
        array = (self.Model - self.Model.mean())*(self.Ref - self.Ref.mean())
        return array.mean()

    def correlation(self):
        return self.covariance()/(np.sqrt(self.Model.std()* self.Ref.std()))

    def bias(self):
        '''
        Bias
        ID 4.3
        Calculates the mean error.
        Result of zero does not necessarily indicate low error due to cancellation.
        '''
        return self.diff().mean();

    def MSE(self):
        '''Mean Square Error
        ID 4.4
        Calculates a mean error (in data units squared), which is not effected by cancellation.
        Squaring the data may cause bias towards large events.
        '''
        return (self.diff()**2).mean()

    def RMSE(self):
        ''' Root mean Square Error
        ID 4.5
        MSE error (4.4) except result is returned in the same units as model,
        which is useful for interpretation.
        '''
        return np.sqrt(self.MSE())

    def MAE(self):
        ''' Mean Absolute Error
        ID 4.6
        Similar to RMSE (4.5) except absolute value is used instead. This reduces
        the bias towards large events; however, it also produces a non-smooth operator
        when used in optimisation.
        '''
        return (np.abs(self.diff())).mean()

    def AME(self):
        '''Absolute Maximum Error
        ID 4.7
        Records the maximum absolute error.
        '''
        return (np.abs(self.diff())).max()



    def NSE(self):
        '''Coefficient of determination
        Nash - Sutcliff Model Efficiency
        ID 6.1

        This method compares the performance of the model to a model that
        only uses the mean of the observed data.
        A value of 1 would indicate a perfect model, while a value of zero
        indicates performance no better than simply using the mean.
        A negative value indicates even worse performance.
        '''

        return 1.0 - self.MSE()/(self.Ref.std()**2)

    def PPMC(self):
        '''
        ID 6.3
        The Pearson Product moment correlation measures the correlation
        of the measured and modelled values.
        Negatives to this model are linear model assumptions and the fact
        it can return an ideal result for a model with constant offset.
        '''
        return 1

    def RSqr(self):
        '''
        Coefficient of determination
        ID 6.4
        Squared version of 6.3, with the same interpretation of results, except range is now (0, 1).
        '''
        return 1

    def IoA(self):
        '''
        Index Of Agreement
        ID 6.5
        This method compares the sum of squared error to the potential error.
        This method is similar to 6.4 however it is designed to be better at handling
        differences in modelled and observed means and variances.
        Squared differences may add bias to large data value events (Willmott, 1981).
        '''

        array_denom = np.abs(self.Model - self.Ref.mean()) + np.abs(self.Model + self.Ref.mean())
        denom = (array_denom**2).sum()
        return 1.0 - self.number()*self.MSE()/denom

    def PI(self):
        '''Persistence Index
        ID 6.6
        The persistence index compares the sum of squared error to the error
        that would occur if the value was forecast as the previous observed value.
        Similar to 6.1 except the performance of the model is being compared to the previous value.
        '''
        s = 0
        for i in range(2,self.number()):
            s+=(self.Ref[i] - self.Ref[i-1])**2
        norm = s/self.number()
        return 1.0 - self.MSE()/norm


    def RAE(self):
        '''Relative Absolute Error
        This compares the total error relative to what the total error
        would be if the mean was used for the model.
        A lower value indicates a better performance, while a score greater
        than one indicates the model is outperformed by using the mean as the prediction.
        '''
        norm = np.abs(self.Ref - self.Ref.mean()).mean()
        return self.MAE()/norm

    def RSR(self):
        ''' RMSE - Standard Deviation  ratio
        The traditional RMSE method weighted by the standard deviation
        of the observed values (Moriasi et al., 2007)
        '''
        return self.RMSE()/self.Ref.std()


class ProfilesMatchup(matchup):
    def __init__(self, Model=None, Ref=None, Depth=None, Lon=None, Lat=None):        
        if Model is None:
            self.Model = np.array([],np.float32)
            self.Ref   = np.array([],np.float32)
            self.Depth = np.array([],np.float32)
            self.Lon   = np.array([],np.float32)
            self.Lat   = np.array([],np.float32)
            self.Lengths = []
        else:
            self.Model = Model
            self.Ref   = Ref
            self.Depth = Depth
            self.Lon   = Lon
            self.Lat   = Lat
            self.Lengths = [len(Model)]
        
    def subset(self,layer):
        if self.number() ==0:
            return self
        ii = (self.Depth <= layer.bottom) & (self.Depth >= layer.top)        
        return FloatMatchup(self.Model[ii], self.Ref[ii], self.Depth[ii], self.Lon[ii], self.Lat[ii])

    
    def extend(self,fm):
        self.Model = np.concatenate((self.Model, fm.Model))
        self.Ref   = np.concatenate((self.Ref,   fm.Ref))
        self.Depth = np.concatenate((self.Depth, fm.Depth))
        self.Lon   = np.concatenate((self.Lon,   fm.Lon))
        self.Lat   = np.concatenate((self.Lat,   fm.Lat))
        self.Lengths.extend([ len(fm.Model)])
    def plot(self):
        '''
        Red lines are reference (biofloat)
        Blue lines are model
        '''

        pl.figure()
        StartInd=0
        for length in self.Lengths:
            End_Ind = StartInd + length
            model = self.Model[StartInd:End_Ind]
            ref   = self.Ref[StartInd:End_Ind]
            pres  = self.Depth[StartInd:End_Ind]

            pl.plot(model,pres,'b', ref,pres,'r')
            StartInd = End_Ind
        pl.gca().invert_yaxis()
        pl.show(block=False)


class ProfileMatchup():
    def __init__(self, Model, Ref, Depth, profileObj):
        bads = np.isnan(Model)
        if bads.any() :
            print "Nans in model "
            Model = Model[~bads]
            Ref   = Ref  [~bads]
            Depth = Depth[~bads]

        self.Model = Model
        self.Ref   = Ref
        self.Depth = Depth
        self.instrument = profileObj
        self.Lon   = np.ones_like(Model)*profileObj.lon
        self.Lat   = np.ones_like(Model)*profileObj.lat

    def plot(self):
        '''
        Red line is reference (biofloat, mooring or vessel)
        Blue line is model
        '''
        pl.figure()
        pl.plot(self.Model,self.Depth,'b', self.Ref,self.Depth,'r')
        pl.gca().invert_yaxis()
        pl.show(block=False)

            

if __name__ == '__main__':
    import numpy as np
    Depth = np.arange(10.)/10
    a=FloatMatchup(np.arange(10),np.arange(10)+1,Depth,Depth,Depth  )
    L=Layer(0.3, 0.8)
    b = a.subset(L)

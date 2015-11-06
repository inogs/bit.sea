import numpy as np
import pylab as pl
from commons.layer import Layer
from statistics import matchup



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
        return ProfilesMatchup(self.Model[ii], self.Ref[ii], self.Depth[ii], self.Lon[ii], self.Lat[ii])

    
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
    a=ProfilesMatchup(np.arange(10),np.arange(10)+1,Depth,Depth,Depth  )
    L=Layer(0.3, 0.8)
    b = a.subset(L)

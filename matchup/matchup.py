import numpy as np
import pylab as pl
from commons.layer import Layer
from statistics import matchup


class FloatProfilesMatchup(ProfilesMatchup):
    def __init__(self, Model=None, Ref=None, Depth=None, Lon=None, Lat=None, time=None, qc=None, name=None, cycle=None):        
        if Model is None:
            self.Model = np.array([],np.float32)
            self.Ref   = np.array([],np.float32)
            self.Depth = np.array([],np.float32)
            self.Lon   = np.array([],np.float32)
            self.Lat   = np.array([],np.float32)
            self.Time  = np.array([],np.float32)
            self.Qc    = np.array([],np.float32)
            self.name  = np.array([],np.float32)
            self.cycle = np.array([],np.float32)
            self.Lengths = []
        else:
            self.Model = Model
            self.Ref   = Ref
            self.Depth = Depth
            self.Lon   = Lon
            self.Lat   = Lat
            self.Time  = time
            self.Qc    = qc
            self.name  = name
            self.cycle = cycle
            self.Lengths = [len(Model)]


    def extend(self,fm):
        self.Model = np.concatenate((self.Model, fm.Model))
        self.Ref   = np.concatenate((self.Ref,   fm.Ref))
        self.Depth = np.concatenate((self.Depth, fm.Depth))
        lenfm = fm.Model.size
        self.Lon   = np.concatenate((self.Lon,   np.ones(lenfm,)*fm.instrument.lon))
        self.Lat   = np.concatenate((self.Lat,   np.ones(lenfm,)*fm.instrument.lat))
        t = fm.instrument.time
        instrumenttime = t.toordinal()+366 + float(t.hour)/24 + float(t.minute)/(24.*60)

        self.Time  = np.concatenate((self.Time,  np.ones(lenfm,)*instrumenttime))
        self.Qc    = np.concatenate((self.Qc,    fm.Qc ))
        self.name  = np.concatenate((self.name,  np.ones(lenfm,)*int(fm.instrument.name()) ))
        self.cycle = np.concatenate((self.cycle, np.ones(lenfm,)*fm.instrument._my_float.cycle))
        self.Lengths.extend([ len(fm.Model)])

    def export(self,directory,prefix):

        self.Model.tofile(directory + prefix + "model.txt",sep="\n",format="%10.5f")
        self.Ref.tofile(  directory + prefix + "ref.txt"  ,sep="\n",format="%10.5f")
        self.Lon.tofile(  directory + prefix + "lon.txt"  ,sep="\n",format="%10.5f")
        self.Lat.tofile(  directory + prefix + "lat.txt"  ,sep="\n",format="%10.5f")
        self.Depth.tofile(directory + prefix + "depth.txt",sep="\n",format="%10.5f")
        self.Time.tofile( directory + prefix + "time.txt" ,sep="\n",format="%10.5f")
        self.Qc.tofile(   directory + prefix + "Qc.txt"   ,sep="\n",format="%10.5f")
        self.name.tofile( directory + prefix + "name.txt" ,sep="\n",format="%d")
        self.cycle.tofile(directory + prefix + "cycle.txt",sep="\n",format="%d")


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
    def __init__(self, Model, Ref, Depth, Qc, profileObj):
        bads = np.isnan(Model)
        if bads.any() :
            print "Nans in model "
            Model = Model[~bads]
            Ref   = Ref  [~bads]
            Depth = Depth[~bads]

        self.Model = Model
        self.Ref   = Ref
        self.Depth = Depth
        self.Qc    = Qc
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

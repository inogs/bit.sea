import numpy as np
import matplotlib.pyplot as pl
from commons.layer import Layer
from matchup.statistics import matchup




class ProfilesMatchup(matchup):
#    def __init__(self, Model=None, Ref=None, Depth=None, Lon=None, Lat=None, Time=None, Qc=None):
    def __init__(self, Model=None, Ref=None, Depth=None, Lon=None, Lat=None, Time=None, Qc=None, name=None):
        if Model is None:
            self.Model = np.array([],np.float32)
            self.Ref   = np.array([],np.float32)
            self.Depth = np.array([],np.float32)
            self.Lon   = np.array([],np.float32)
            self.Lat   = np.array([],np.float32)
            self.Time  = np.array([],np.float32)
            self.Qc    = np.array([],np.float32)
            self.name  = []
            self.Lengths = []
            self.CheckReports=[]
        else:
            self.Model = Model
            self.Ref   = Ref
            self.Depth = Depth
            self.Lon   = Lon
            self.Lat   = Lat
            self.Time  = Time
            self.Qc    = Qc
            self.name  = name
            self.Lengths = [len(Model)]

    def subset(self,layer):
        '''
        Subset using a layer

        Argument:
        *  layer * Layer object

        Returns:
        * P * ProfilesMatchup object, where is  layer.top <= P.Depth <= layer.bottom
        '''
        if self.number() ==0:
            return self
        ii = (self.Depth <= layer.bottom) & (self.Depth >= layer.top)
        names = []
        for iin,iname in enumerate(ii): 
            if iname: names.append(self.name[iin])
        return ProfilesMatchup(self.Model[ii], self.Ref[ii], self.Depth[ii], self.Lon[ii], self.Lat[ii], self.Time[ii], self.Qc[ii], names)

    def limmaxref(self,value):
        '''
        Subset about Ref data

        Argument:
        *  value * float, maximum value for Ref

        Returns:
        * P * ProfilesMatchup object, where is P.Ref <= value
        '''
        if self.number() ==0:
            return self
        ii = (self.Ref <= value)
        names = []
        for iin,iname in enumerate(ii): 
            if iname: names.append(self.name[iin])
        return ProfilesMatchup(self.Model[ii], self.Ref[ii], self.Depth[ii], self.Lon[ii], self.Lat[ii], self.Time[ii], self.Qc[ii], names)

    def limminmodel(self,value):
        '''
        Subset about Model data

        Argument:
        *  value * float, minimum value for Model

        Returns:
        * P * ProfilesMatchup object, where is P.Model >= value
        '''
        if self.number() ==0:
            return self
        ii = (self.Model >= value)
        names = []
        for iin,iname in enumerate(ii): 
            if iname: names.append(self.name[iin])
        return ProfilesMatchup(self.Model[ii], self.Ref[ii], self.Depth[ii], self.Lon[ii], self.Lat[ii], self.Time[ii], self.Qc[ii], names)

    def extend(self,fm):
        '''
        Extends itself without returning anything.

        Argument:
        * fm * a ProfileMatchup or ProfilseMatchup object

        '''
        lenfm = fm.Model.size
        fmLon,fmLat,fmTime = fm.localizationArrays()

        self.Model = np.concatenate((self.Model, fm.Model))
        self.Ref   = np.concatenate((self.Ref,   fm.Ref))
        self.Depth = np.concatenate((self.Depth, fm.Depth))
        self.Lon   = np.concatenate((self.Lon,   fmLon))
        self.Lat   = np.concatenate((self.Lat,   fmLat))
        self.Time  = np.concatenate((self.Time,  fmTime))
        self.Qc    = np.concatenate((self.Qc,    fm.Qc ))
        self.name.extend( [fm.instrument.name() for k in range(lenfm)] )

        self.Lengths.extend([ lenfm ])
        self.CheckReports.append(fm.checkreport)

    def export(self,directory,prefix):

        self.Model.tofile(directory + prefix + "model.txt",sep="\n",format="%10.5f")
        self.Ref.tofile(  directory + prefix + "ref.txt"  ,sep="\n",format="%10.5f")
        self.Lon.tofile(  directory + prefix + "lon.txt"  ,sep="\n",format="%10.5f")
        self.Lat.tofile(  directory + prefix + "lat.txt"  ,sep="\n",format="%10.5f")
        self.Depth.tofile(directory + prefix + "depth.txt",sep="\n",format="%10.5f")
        self.Time.tofile( directory + prefix + "time.txt" ,sep="\n",format="%10.5f")
        self.Qc.tofile(   directory + prefix + "Qc.txt"   ,sep="\n",format="%10.5f")
        np.array(self.name).tofile(directory + prefix + "name.txt" ,sep="\n",format="%s")


    def plot(self,fig=None,ax=None):
        '''
        Red lines are reference (biofloat)
        Blue lines are model
        '''
        if (fig is None) or (ax is None):
            fig , ax = pl.subplots()
        StartInd=0
        for il, length in enumerate(self.Lengths):
            End_Ind = StartInd + length
            model = self.Model[StartInd:End_Ind]
            ref   = self.Ref[StartInd:End_Ind]
            pres  = self.Depth[StartInd:End_Ind]
            if il==0:
                ax.plot(model,pres,'b', label='model')
                ax.plot(ref,pres,'r', label='float')
            else:
                ax.plot(model,pres,'b',ref,pres,'r')
            StartInd = End_Ind
        if not ax.yaxis_inverted(): ax.invert_yaxis()
        ax.legend()
        return fig,ax

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
        lenfm = fm.Model.size
        fmLon,fmLat,fmTime = fm.localizationArrays()

        self.Model = np.concatenate((self.Model, fm.Model))
        self.Ref   = np.concatenate((self.Ref,   fm.Ref))
        self.Depth = np.concatenate((self.Depth, fm.Depth))
        self.Lon   = np.concatenate((self.Lon,   fmLon))
        self.Lat   = np.concatenate((self.Lat,   fmLat))
        self.Time  = np.concatenate((self.Time,  fmTime))
        self.Qc    = np.concatenate((self.Qc,    fm.Qc ))

        self.name  = np.concatenate((self.name,  np.ones(lenfm,)*int(fm.instrument.name()) ))
        self.cycle = np.concatenate((self.cycle, np.ones(lenfm,)*fm.instrument._my_float.cycle))
        self.Lengths.extend([ lenfm ])

    def export(self,directory,prefix):

        self.Model.tofile(directory + prefix + "model.txt",sep="\n",format="%10.5f")
        self.Ref.tofile(  directory + prefix + "ref.txt"  ,sep="\n",format="%10.5f")
        self.Lon.tofile(  directory + prefix + "lon.txt"  ,sep="\n",format="%10.5f")
        self.Lat.tofile(  directory + prefix + "lat.txt"  ,sep="\n",format="%10.5f")
        self.Depth.tofile(directory + prefix + "depth.txt",sep="\n",format="%10.5f")
        self.Time.tofile( directory + prefix + "time.txt" ,sep="\n",format="%10.5f")
        self.Qc.tofile(   directory + prefix + "Qc.txt"   ,sep="\n",format="%10.5f")
        np.array(self.name).tofile(directory + prefix + "name.txt" ,sep="\n",format="%s")
        self.cycle.tofile(directory + prefix + "cycle.txt",sep="\n",format="%d")


class ProfileMatchup():
    def __init__(self, Model, Ref, Depth, Qc, profileObj, checkreport=None, accept_nans=False):
        bads = np.isnan(Model)
        if not accept_nans:
            if bads.any() :
                print("matchup: Nans in model ")
                Model = Model[~bads]
                Ref   = Ref  [~bads]
                Depth = Depth[~bads]
                Qc    = Qc[~bads]

        self.Model = Model
        self.Ref   = Ref
        self.Depth = Depth
        self.Qc    = Qc
        self.instrument = profileObj
        self.checkreport =checkreport


    def localizationArrays(self):
        '''Returns arrays Lon,Lat,Time as repetition of instrument fields
        Time is expressed as matlab time
        '''
        n = len(self.Model)
        Lon   = np.ones(n,)*self.instrument.lon
        Lat   = np.ones(n,)*self.instrument.lat
        t = self.instrument.time
        instrumenttime = t.toordinal()+366 + float(t.hour)/24 + float(t.minute)/(24.*60)
        Time  = np.ones(n,)*instrumenttime
        return Lon,Lat,Time

    def plot(self):
        '''
        Red line is reference (biofloat, mooring or vessel)
        Blue line is model
        '''
        pl.figure()
        pl.plot(self.Model,self.Depth,'b', self.Ref,self.Depth,'r')
        pl.gca().invert_yaxis()
        pl.show(block=False)


    def plot_file(self,date,p,var,dirout):
        '''
        Red line is reference (biofloat, mooring or vessel)
        Blue line is model
        '''
        floatname = p.name()
        filename =dirout+"/"+date.strftime('%Y%m%d') +"_"+floatname+"_"+var+ ".png"
        pl.figure()
        pl.plot(self.Model,self.Depth,label="Model")
        pl.plot(self.Ref,self.Depth,label="Float")
        pl.legend(loc='upper right')
        pl.ylabel("depth")
        figtitle = var+" date="+date.strftime('%Y/%m/%d')+" float="+floatname
        pl.title(figtitle)
        pl.gca().invert_yaxis()
        pl.savefig(filename)
        pl.close()
        import libxmp
        import  libxmp.utils
        from libxmp import XMPFiles, consts
        xmpfile = XMPFiles( file_path=filename, open_forupdate=True )
        xmp = xmpfile.get_xmp()

        if xmp is None:
            xmp = libxmp.XMPMeta()

        xmp.set_property(consts.XMP_NS_DC, 'float', floatname )
        xmp.set_property(consts.XMP_NS_DC, 'date', date.strftime('%Y%m%d') )
        xmp.set_property(consts.XMP_NS_DC, 'hour', date.strftime('%H:%M:%S') )
        xmp.set_property(consts.XMP_NS_DC, 'var', var )
        xmp.set_property(consts.XMP_NS_DC, 'position.lat',str(p.lat)+"N")
        xmp.set_property(consts.XMP_NS_DC, 'position.lon',str(p.lon)+"E")
        xmpfile.put_xmp(xmp)
        xmpfile.close_file()

    def plot_subplot(self,var,fig, ax):
        '''
        Red line is reference (biofloat, mooring or vessel)
        Blue line is model
        '''
        import latexcodec
        ax.plot(self.Model,self.Depth,'r')
        ax.plot(self.Ref,self.Depth,'b')


        figtitle = var
        ax.set_title(figtitle.encode("latex"))
        ax.invert_yaxis()

        return fig, ax



if __name__ == '__main__':
    Depth = np.arange(10.)/10
    a=ProfilesMatchup(np.arange(10),np.arange(10)+1,Depth,Depth,Depth,Depth,Depth  )
    L=Layer(0.3, 0.8)
    b = a.subset(L)

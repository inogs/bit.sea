import pickle
import matplotlib.pyplot as pl
from basins import V2 as OGS

class filereader():
    def __init__(self, filename):
        
        fid = open(filename,'rb')
        LIST = pickle.load(fid)
        fid.close()
        TIMES,_,_,MODEL_MEAN,SAT___MEAN,_,_,MODEL__STD,SAT____STD,CORR = LIST
        self.TIMES = TIMES
        self.MODEL_MEAN = MODEL_MEAN
        self.SAT___MEAN  = SAT___MEAN
        self.MODEL__STD  = MODEL__STD
        self.SAT____STD  = SAT____STD
        


INPUTDIR="/Users/gbolzon/Documents/workspace/bit.sea/"

VARLIST= ["P1l","P2l","P3l","P4l"]

MATRIX_LIST=[filereader(INPUTDIR + var  +'_open_sea.pkl') for var in VARLIST]

#P1l = filereader(INPUTDIR + 'P1l_open_sea.pkl')
#P2l = filereader(INPUTDIR + 'P2l_open_sea.pkl')
#P3l = filereader(INPUTDIR + 'P3l_open_sea.pkl')
#P4l = filereader(INPUTDIR + 'P4l_open_sea.pkl')



isub = 0
sub = OGS.alb
pl.close('all')
fig = pl.figure(figsize=(10, 10))
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)

AXES_LIST=[ax1, ax2, ax3, ax4]


model_label='MODEL'
units="[mg/m$^3$]"
var_label = 'Sat' + " " + units


for ivar,var in enumerate(VARLIST):
    ax = AXES_LIST[ivar]    
    Pl = MATRIX_LIST[ivar]
    ax.plot(Pl.TIMES,Pl.SAT___MEAN[:,isub],'og',label=' SAT')
    ax.fill_between(Pl.TIMES,Pl.SAT___MEAN[:,isub]-Pl.SAT____STD[:,isub],Pl.SAT___MEAN[:,isub]+Pl.SAT____STD[:,isub],color='palegreen')
    ax.plot(Pl.TIMES,Pl.MODEL_MEAN[:,isub],'-k',label=model_label)
    ax.plot(Pl.TIMES,Pl.MODEL_MEAN[:,isub]-Pl.MODEL__STD[:,isub],':k')
    ax.plot(Pl.TIMES,Pl.MODEL_MEAN[:,isub]+Pl.MODEL__STD[:,isub],':k')
    ax.grid(True)
    ax.set_ylabel("%s - %s" %(var,units))
    if ax != ax4: ax.set_xticklabels([])

ax1.set_title("%s" %(sub.name.upper()  ) ).set_fontsize(14)
ax1.legend()
xlabels = ax4.get_xticklabels()
pl.setp(xlabels, rotation=30)    

fig.show()
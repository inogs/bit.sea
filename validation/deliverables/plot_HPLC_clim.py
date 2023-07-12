import argparse
from commons.utils import addsep
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates plots comparison of 4PFTs with HPLC dataset
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--inputdir','-i',
                                type = str,
                                required = True,
                                help = '')
    parser.add_argument(   '--outdir','-o',
                                type = str,
                                required = True,
                                help = '')

    parser.add_argument(   '--ptype','-p',
                                type = str,
                                required = True,
                                help = 'Phytoplankton Functional Type')
    return parser.parse_args()

import matplotlib
matplotlib.use('Agg')

import pandas as pd
from datetime import datetime
import csv
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

args = argument()

INDIR=addsep(args.inputdir)
OUTDIR=addsep(args.outdir)
P_type=args.ptype

#INDIR="/g100_scratch/userexternal/lfeudale/validation/V10C/run4.19/bit.sea/validation/deliverables/"
#P_type="P4l"

if (P_type=="P1l"): P_name="DIATO"
if (P_type=="P2l"): P_name="NANO"
if (P_type=="P3l"): P_name="PICO"
if (P_type=="P4l"): P_name="DINO"

infile_MEAN=INDIR + P_type + "_OPTION3_newSeasons3_means.csv"
infile_STD=INDIR  + P_type + "_OPTION3_newSeasons3_std.csv"

headers = ["LAYER","meanWref","meanWmod","meanSref","meanSmod"]
headers_STD = ["LAYER","stdWref","stdWmod","stdSref","stdSmod"]
layers = ["0-10","10-30","30-60","60-100 ","100-150","150-300","0-300"]
layer_mean_depth = ["5","20","50","80","125","225","150"]
depths = [5,20,50,80,125,225,150]

#Pl1=pd.read_csv(fileP1,sep=",",index_col=0)   
Pl=pd.read_csv(infile_MEAN,sep=",",index_col=0,usecols=headers) 
Pl_std=pd.read_csv(infile_STD,sep=",",index_col=0,usecols=headers_STD)
#Pl1=pd.read_csv(fileP1,delim_whitespace = True, index_col=0, usecols=headers)

fig,(ax1,ax2) =plt.subplots(1,2) 
ax1.plot(Pl.iloc[:-1,0],layers[:-1],':',color="blueviolet" ) 
ax1.plot(Pl.iloc[:-1,1],layers[:-1],'-',color="cyan" )  
ax2.plot(Pl.iloc[:-1,2],layers[:-1],':',color="coral" ) 
ax2.plot(Pl.iloc[:-1,3],layers[:-1],'-',color="orange" ) 

#ax.fill_between(Pl.iloc[:-1,0]-Pl_std[:-1,0],Pl.iloc[:-1,0]+Pl_std[:-1,0],layers[:-1],color=lavender)
#ax.fill_between(layers[:-1],Pl.iloc[:-1,0]-Pl_std.iloc[:-1,0],Pl.iloc[:-1,0]+Pl_std.iloc[:-1,0],color="lavender")
ax1.fill_betweenx(layers[:-1],Pl.iloc[:-1,0]-Pl_std.iloc[:-1,0],Pl.iloc[:-1,0]+Pl_std.iloc[:-1,0],color="lavender")
ax1.fill_betweenx(layers[:-1],Pl.iloc[:-1,1]-Pl_std.iloc[:-1,1],Pl.iloc[:-1,1]+Pl_std.iloc[:-1,1],color="lightcyan")
ax2.fill_betweenx(layers[:-1],Pl.iloc[:-1,2]-Pl_std.iloc[:-1,2],Pl.iloc[:-1,2]+Pl_std.iloc[:-1,2],color="mistyrose")
ax2.fill_betweenx(layers[:-1],Pl.iloc[:-1,3]-Pl_std.iloc[:-1,3],Pl.iloc[:-1,3]+Pl_std.iloc[:-1,3],color="moccasin")

ax1.set_xlim([0,0.25])
ax2.set_xlim([0,0.25])
#ax2.set(ylabel=None)

# Set the tick positions
ax1.set_yticks(layers[:-1])
# Set the tick labels
ax1.set_yticklabels(layer_mean_depth[:-1])
ax2.set_yticklabels([])
#ax1.set_xlabel(P_type + " " + " $[ mg {\,} m^{-3} ]$")
#ax2.set_xlabel(P_type + " " + " $[ mg {\,} m^{-3} ]$")
ax1.set_xlabel("$[ mg {\,} m^{-3} ]$")
ax2.set_xlabel("$[ mg {\,} m^{-3} ]$")

ax1.set_ylabel("depth [m]")


ax1.invert_yaxis() 
ax2.invert_yaxis()
#ax.legend(['WIN ref','WIN mod','SUM ref','SUM mod'])  
ax1.legend(['HPLC','mod'])
ax2.legend(['HPLC','mod'])
ax1.set_title("WINTER")
ax2.set_title("SUMMER")
fig.suptitle(P_name)
#ax.invert_xaxis()                                                                                               
#ax.invert_yaxis() 
fig.savefig(OUTDIR + 'HPLC_profiles_comparison_STD_' + P_type + '_seasons.png')
#plt.savefig('HPLC_profiles.png') 

plt.close(fig) 

#for ilay, lay in enumerate(layers):
#   df1=pd.DataFrame({'Wr':P1l.iloc[ilay,0],'Wm':P1l.iloc[ilay,1],'Sr':P1l.iloc[ilay,2],'Sm':P1l.iloc[ilay,3]},index=headers )

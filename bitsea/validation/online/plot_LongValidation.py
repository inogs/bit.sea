import argparse

def argument():
    parser = argparse.ArgumentParser(description='''
    plot validation time series
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--indir', '-i',
                        type=str,
                        default=None,
                        required=True,
                        help="Dir with hc validation")
    
    parser.add_argument('--outdir', '-o',
                        type=str,
                        default=None,
                        required=True,
                        help="output dir")

    return parser.parse_args()

args = argument()


import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as md
from commons.utils import addsep



INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)


def extractLIST(infile):

    table = np.loadtxt(infile,dtype='str')
    ndates = table.shape[0]

    LISTdates = [] 
    LISTrmsd = []

    for ii in range(ndates):
        dd = table[ii,0]
        date = datetime.datetime.strptime(dd,'%Y-%m-%d')
        LISTdates.append(date)
        LISTrmsd.append(float(table[ii,2]))

    return LISTdates,LISTrmsd


infile = INDIR + 'table_statistics_V9C_DT.txt'
LISTdates, LISTrmsd = extractLIST(infile)


LISTfilermsd_old = [
    'table_statistics_MERSEA.txt',
    'table_statistics_MyOcean1.txt',
    'table_statistics_MyOcean2.txt',
    'table_statistics_Copernicus_till_V8C_DT.txt'
]

LISTdates_old = {}
LISTrmsd_old = {}

for pp in LISTfilermsd_old:
    infile_old = INDIR + '/' + pp
    LISTdates_old[pp],LISTrmsd_old[pp] = extractLIST(infile_old)



LISTarrows = [
    'table_arrow_MERSEA.txt',
    'table_arrow_MyOcean.txt',
    'table_arrow_Copernicus.txt',
]

DICTarrnames = {
    'table_arrow_MERSEA.txt': 'MERSEA',
    'table_arrow_MyOcean.txt': 'MyOcean',
    'table_arrow_Copernicus.txt': 'Copernicus Marine Service',
}

LISTdates_arr = {}
LISTy_arr = {}

for pp in LISTarrows:
    infile_arr = INDIR + '/' + pp
    LISTdates_arr[pp],LISTy_arr[pp] = extractLIST(infile_arr)



LISTresol = [
    'table_res8.txt',
    'table_res16.txt',
    'table_res24.txt',
]

DICTresnames = {
    'table_res8.txt': '1/8°',
    'table_res16.txt': '1/16°',
    'table_res24.txt': '1/24°',
}

LISTdates_res = {}
LISTy_res = {}

for pp in LISTresol:
    infile_res = INDIR + '/' + pp
    LISTdates_res[pp],LISTy_res[pp] = extractLIST(infile_res)


LISTevents = [
    'table_event_DA.txt',
    'table_event_BC.txt',
    'table_event_Dcycle.txt',
    'table_event_Optics.txt',
]

Nevents = len(LISTevents)

DICTevnames = {
    'table_event_DA.txt': 'BGC DA',
    'table_event_BC.txt': 'New Med-PHY BC at Gibraltar',
    'table_event_Dcycle.txt': 'Daily cycle in BFM',
    'table_event_Optics.txt': 'Optics',
}

DICTevcolors = {
    'table_event_DA.txt': 'purple',
    'table_event_BC.txt': 'purple',
    'table_event_Dcycle.txt': 'purple',
}
DICTevcolors = {
    'table_event_DA.txt': 'orchid',
    'table_event_BC.txt': 'mediumpurple',
    'table_event_Dcycle.txt': 'royalblue',
    'table_event_Optics.txt': 'teal',
}

LISTdates_ev = {}
LISTy_ev = {}

for pp in LISTevents:
    infile_ev = INDIR + '/' + pp
    LISTdates_ev[pp],LISTy_ev[pp] = extractLIST(infile_ev)



plt.close('all')

fig,ax = plt.subplots(1,1,figsize=[10,5])


for pp in LISTfilermsd_old:
    plt.plot(LISTdates_old[pp],LISTrmsd_old[pp],'--',color='salmon')
    if ('Copernicus' in pp)==False:
        plt.plot(LISTdates_old[pp][0],LISTrmsd_old[pp][0],'o',color='r',
            markersize=12,
            alpha=0.3
            )
        plt.plot(LISTdates_old[pp][0],LISTrmsd_old[pp][0],'o',color='r',
            markersize=12,
            mfc='None'
            )

plt.plot(LISTdates,LISTrmsd,'o-',color='r',
        markersize=3,
        linewidth=.3,
        )

xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
# ax.set_ylim(0,ymax+ymax*.1)


for pp in LISTarrows:
    datenum0 = md.date2num(LISTdates_arr[pp][0])
    datenum1 = md.date2num(LISTdates_arr[pp][1])
    if datenum1>xmax: datenum1=xmax

    plt.annotate(DICTarrnames[pp],xy=(datenum1,LISTy_arr[pp][1]),
        xytext=(datenum0,LISTy_arr[pp][0]),
        color='orange',
        horizontalalignment='left',
        verticalalignment='bottom',
        arrowprops=dict(arrowstyle="-|>",color='orange',linewidth='3',relpos=(-0.1,0))
        )

for pp in LISTresol:
    datenum0 = md.date2num(LISTdates_res[pp][0])
    datenum1 = md.date2num(LISTdates_res[pp][1])
    if datenum1>xmax: datenum1=xmax
    if LISTy_res[pp][1]<ymax:
        LISTy_res[pp][1] = ymax
    if LISTy_res[pp][0]<ymax:
        LISTy_res[pp][0] = ymax

    plt.axvspan(datenum0-60,datenum0,color='teal',alpha=0.15)

    plt.annotate(DICTresnames[pp],xy=(datenum1,LISTy_res[pp][1]),
        xytext=(datenum0,LISTy_res[pp][0]),
        color='teal',
        horizontalalignment='left',
        verticalalignment='bottom',
        arrowprops=dict(arrowstyle="-|>",color='teal',linewidth='3',relpos=(-0.1,0))
        )


delta = 0.05
shift = 1.1

for ii,pp in enumerate(LISTevents):
    datenum0 = md.date2num(LISTdates_ev[pp][0])
    datenum1 = md.date2num(LISTdates_ev[pp][1])
    if datenum1>xmax: datenum1=xmax

    xratio = (datenum0-xmin)/(xmax-xmin)
    if ii%2==0:
        val = 'bottom'
        idelta = 1.1
    else:
        val = 'top'
        idelta = -0.2
    plt.axhspan(ymax*(shift),ymax*(shift+delta),xmin=xratio,
            # color=DICTevcolors[pp],alpha=0.2)
            color=DICTevcolors[pp])
    plt.text(datenum0,ymax*(shift+idelta*delta),DICTevnames[pp],
            va=val,ha='left',color=DICTevcolors[pp])




xmin = md.date2num(datetime.datetime(2007,3,1))
ax.set_xlim(xmin,xmax)
# ax.set_ylim(0,ymax+delta*(Nevents))
ax.set_ylim(0,ymax+delta*2.5)

plt.ylabel('Log' + r'$_{10}$' + ' [mg chl/m' + r'$^3$' + ']')

plt.grid()

plt.tight_layout()

#plt.show(block=False)


plt.savefig(OUTDIR + 'longRMSD.png')



import pickle
import matplotlib.pyplot as pl
import matplotlib.dates as mdates

fid = open('export_data_ScMYValidation_plan.pkl')
LIST = pickle.load(fid)
fid.close()
TIMES                          = LIST[0]
BGC_CLASS4_CHL_RMS_SURF_BASIN  = LIST[1]
BGC_CLASS4_CHL_BIAS_SURF_BASIN = LIST[2]
MODEL_MEAN                     = LIST[3]
SAT___MEAN                     = LIST[4]



from basins import OGS 
for isub,sub in enumerate(OGS.P): 
    print sub.name
    fig, ax = pl.subplots()
    fig.set_dpi = 100
    ax.plot(TIMES,MODEL_MEAN[:,isub])
    ax.set_ylabel('CHL')
    ax.xaxis.set_major_locator(mdates.YearLocator())
    #fig.show()
    outfilename = "chl" + sub.name + ".png"
    #pl.savefig()
    import sys
    sys.exit()



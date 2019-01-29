from basins import V2 as OGS
from static import socat_reader
import numpy as np
from commons import timerequestors
from commons.utils import writetable
S=socat_reader.CO2_socat_reader()


nSUB = len(OGS.P.basin_list)
CLIM=np.zeros((nSUB, 12),np.float32)*np.nan
CLIM_STD=np.zeros((nSUB, 12),np.float32)*np.nan
nP_month=np.zeros((nSUB, 12),np.int32)*np.nan

for isub, sub in enumerate(OGS.P):
    print sub
    for imonth in range(12):
        print imonth
        monthly_req= timerequestors.Clim_month(imonth+1)
        values=S.clim_month_selector('pCO2', imonth+1, sub)
        nP = len(values)
        if nP >0:
            V=np.array(values)
            CLIM[isub,imonth]=V.mean()
            CLIM_STD[isub,imonth]=V.std()
            nP_month[isub,imonth]=nP
            print nP
#         Profilelist=S.Selector('fCO2', monthly_req, sub)
#         nP = len(Profilelist)
#         if nP >0:
#             values=np.zeros((nP,),np.float32)
#             for ip, p in enumerate(Profilelist):
#                 Pres,profile,Qc=p.read('fCO2')
#                 values[ip]=profile.mean() # ci puo essere piu di un valore
#             
#             CLIM[isub,imonth]=values.mean()

rows_names_list=[sub.name for sub in OGS.P]
column_names_list=[str(i) for i in range(1,13)]
writetable("monthly_clim_socat.txt", CLIM, rows_names_list, column_names_list)
writetable("monthly_clim_socat_STD.txt", CLIM_STD, rows_names_list, column_names_list)
writetable("monthly_num_socat.txt", nP_month, rows_names_list, column_names_list)

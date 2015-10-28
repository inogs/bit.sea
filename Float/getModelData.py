from shared_data import *
import MatchGenerator
from commons.time_interval import TimeInterval

T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d-%H:%M:%S')

M = MatchGenerator.Float_Matchup_Manager(T_INT,INPUTDIR,BASEDIR)
M.writefiles_for_profiling('./jobProfiler.sh')
M.dumpModelProfiles('./jobProfiler.sh')





# Coupled_List=[]
# for ir, req in enumerate(TL.getOwnList()):
#     LIST_of_REQ=[]
#     for f in FLOAT_LIST:
#         if (f.time >= req.starttime) & (f.time <= req.endtime) :
#             LIST_of_REQ.append(f)
#     if (len(LIST_of_REQ) >0 ): Coupled_List.append((TL.Timelist[ir],LIST_of_REQ))




from shared_data import *
import MatchGenerator

M = MatchGenerator.Float_Matchup_Manager(DATESTART,DATE__END,INPUTDIR)
M.writefiles_for_profiling('jobProfiler.sh',BASEDIR)
M.dumpModelProfiles('jobProfiler.sh')





# Coupled_List=[]
# for ir, req in enumerate(TL.getOwnList()):
#     LIST_of_REQ=[]
#     for f in FLOAT_LIST:
#         if (f.time >= req.starttime) & (f.time <= req.endtime) :
#             LIST_of_REQ.append(f)
#     if (len(LIST_of_REQ) >0 ): Coupled_List.append((TL.Timelist[ir],LIST_of_REQ))




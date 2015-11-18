
from instruments.matchup_manager import Matchup_Manager
from commons.time_interval import TimeInterval

INPUTDIR="/pico/scratch/userexternal/gbolzon0/PROFILATORE/AVE/"
BASEDIR='/pico/scratch/userexternal/gbolzon0/PROFILATORE/'

DATESTART = '20150901-00:00:00'
DATE__END = '20150917-00:00:00'

T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d-%H:%M:%S')


if __name__ == '__main__':

    M = Matchup_Manager(T_INT,INPUTDIR,BASEDIR)
    M.writefiles_for_profiling('./jobProfiler.sh')
    M.dumpModelProfiles('./jobProfiler.sh')



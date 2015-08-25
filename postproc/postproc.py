import scipy.io.netcdf as NC
import glob
import os, sys
import GB_lib as G
import read_descriptor
import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''Creates big ave files for archive and/or big files with aggregated var for 
    Avescan.py. Creates also chl sup.''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = '/some/path/MODEL/AVE_FREQ_1/')

    parser.add_argument(   '--tmpdir', '-t',
                                type = str,
                                default = None,
                                help = """ /some/path/POSTPROC/output/AVE_FREQ_1/TMP/ .  
                                Path to put files with aggregated variables for aveScan.py. 
                                No generation of aggregate files if this parameter is omitted.
                                """)
    parser.add_argument(   '--archivedir', '-a',
                                type = str,
                                default = None,
                                help = '''/some/path/POSTPROC/output/AVE_FREQ_1/Archive/  . 
                                Path to put native vars as they are, in order to compress them.
                                No generation of archived files if this parameter is omitted.
                                ''')
    
    parser.add_argument(   '--chlsupdir', '-c',
                                type = str,
                                default = None,
                                help = """/some/path/POSTPROC/output/AVE_FREQ_1/CHL_SUP.
                                No generation of chl sup if this parameter is omitted.
                                """)    
    parser.add_argument(   '--avelist',"-l",
                                type = str,
                                default = "ave*N1p.nc",
                                help = 'ave*.N1p.nc, they configure the date list')
    parser.add_argument(   '--descriptor',"-d",
                                type = str,
                                required = True,
                                help = 'VarDescriptor_1.xml, or the complete path')       

    return parser.parse_args()

def addsep(string):    
    if string[-1] != os.sep:
        return string + os.sep
    else:
        return  string
    
        
    

try :
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks = comm.size 
except:
    rank   = 0
    nranks = 1

args = argument()
 

AVEDIR     = addsep(args.inputdir)

RD = read_descriptor.read_descriptor(args.descriptor)
PATH_NAME = AVEDIR + args.avelist
if rank==0 : print "INPUT_DIR =", AVEDIR


if args.archivedir:
    ARCHIVEdir = addsep(args.archivedir)
    if rank==0 : print "ARCHIVEDIR=", ARCHIVEdir
    os.system("mkdir -p " + ARCHIVEdir)

if args.tmpdir:
    TMPOUTdir  = addsep(args.tmpdir)
    if rank==0 : print "TMPOUTDIR= ", TMPOUTdir
    os.system("mkdir -p " + TMPOUTdir)

if args.chlsupdir:
    CHLSUPdir  = addsep(args.chlsupdir)
    if rank==0 : print "CHLSUPDIR =", CHLSUPdir
    os.system("mkdir -p " + CHLSUPdir)

SingleVar_filelist=glob.glob(PATH_NAME)
SingleVar_filelist.sort()

for N1pfile in SingleVar_filelist[rank::nranks]:
    varname   = N1pfile[-6:-3] 
    dailyAve  = os.path.basename(N1pfile).replace("."+ varname,"")    

    print "writing ", dailyAve
    
    if args.archivedir :
        Big_Ave__archive   = ARCHIVEdir + dailyAve
        G.WriteBigAve(N1pfile, Big_Ave__archive, RD.ARCHIVE_VARS)
    
    if args.tmpdir:
        Big_Ave_postproc   = TMPOUTdir + dailyAve
        G.WriteTMPave(N1pfile, Big_Ave_postproc, RD)
    
    if args.chlsupdir:    
        chlfile =dailyAve.replace("ave","chl")
        G.writeChlSup(Big_Ave_postproc, CHLSUPdir+chlfile)

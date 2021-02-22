import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Find differences between two argo floats files
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--input_NEW',"-N",
                                type = str,
                                required = True,
                                help = 'input file argo NEW, eg=')
    parser.add_argument(   '--input_OLD',"-O",
                                type = str,
                                required = True,
                                help = 'input file argo txt OLD')
    parser.add_argument(   '--outputfile',"-o",
                                type = str,
                                required = True,
                                help = 'output file name')

    return parser.parse_args()

args = argument()

from commons.utils import file2stringlist

OLD=file2stringlist(args.input_OLD)
NEW=file2stringlist(args.input_NEW)
outputfile=args.outputfile


DIFF=[]

for line in NEW:
    if line not in OLD:
	DIFF.append(line)


LINES=[]
for line in DIFF:
    LINES.append(line + '\n')

F = file(outputfile,'w')
F.writelines(LINES)
F.close()


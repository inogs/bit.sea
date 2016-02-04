import os
def addsep(string):    
    if string[-1] != os.sep:
        return string + os.sep
    else:
        return  string


def file2stringlist(filename):
    '''
    Argument : string - a path name indicating a text file
    Returns  : a list of strings
  
    A text file is converted in a list of strings, one for each row.
    
    '''
    LIST=[]
    filein=open(filename)
    for line in filein:
        LIST.append(line[:-1])
    filein.close()
    return LIST


def find_index(thestring, STRINGLIST):
    '''
    Searches a string in a list of strings.
    The string is supposed exist in the list. If not, the function raises an error.
    The string list is supposed be like in NetCDF files, where a list of variable names is
    written in an only 2D character array. 
    '''
    nStrings = STRINGLIST.shape[0]
    for istring in range(nStrings):
        strippedstring=STRINGLIST[istring,:].tostring().strip()
        if strippedstring == thestring: break
    else:
        print thestring + " Not Found"
        raise NameError
    return istring
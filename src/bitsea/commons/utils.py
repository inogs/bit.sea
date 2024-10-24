# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import os
import numpy as np
import re
"""Helper functions"""


def get_date_string(s):
    """Finds a date in YYYYMMDD format into a bigger string.

    Args:
        - *s*: the string.

    Returns: a tuple with the date string in YYYY-MM-DD format as first element
    and YYYYMMDD format as second element. If no date has been found returns
    two empty strings.
    """
    m = re.search('.*([0-9]{4})([0-9]{2})([0-9]{2}).*',s)
    if m is None:
        return "", ""
    longdate = "%s-%s-%s" % (m.group(1), m.group(2), m.group(3))
    shortdate = "%s%s%s" % (m.group(1), m.group(2), m.group(3))
    return longdate, shortdate

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
    try:
        return STRINGLIST.index(thestring)
    except IndexError:
        print(thestring + " Not Found")
        raise NameError('Variable should be one of the following: {}'.format(STRINGLIST))

def find_index_s(substring, STRINGLIST):
    '''
    Searches a substring in a list of strings.
    It returns the array of indexes.
    '''
    subst=substring.lower()
    idx=[]
    namesC=[]
    for istring, current_str in enumerate(STRINGLIST):
        if subst in current_str.lower(): 
            # print(current_str.lower())
            idx.append(istring)
            namesC.append(current_str)
    return np.array(idx),namesC


def ticklabels_degree(ax,fsize=7,intdeg=True):
    '''
    Argument : 
         - ax: ax object (e.g. ax = fig.add_subplot())
              the tick labels of ax will be modified
         - fsize: fontsize of the tick labels (default=7)
         - intdeg: if True the label values are integer
                   if False the label values are decimal
                   (default is True)
  
    Modify ticklabels of ax in format degE and degN (deg is the degree symbol).
    '''
    degsymb = u'\xb0'
    xticks = ax.get_xticks()
    xticklabels = []
    for xt in xticks:
        if intdeg: xt = np.int(xt)
        xticklabels.append(np.str(xt) + degsymb + 'E')
    ax.set_xticklabels(xticklabels, fontsize=7)
    yticks = ax.get_yticks()
    yticklabels = []
    for yt in yticks:
        if intdeg: yt = np.int(yt)
        yticklabels.append(np.str(yt) + degsymb + 'N')
    ax.set_yticklabels(yticklabels, fontsize=7)


def getcolor(ntimes,itime, colormap='gist_ncar'):
    '''
    Uses matplotlib colormap
    Arguments:
    * nTimes   * integer, the number of times
    * itime    * integer, time index
    * colormap * string (optional), is the name of matplotlib colormap
                 See https://matplotlib.org/examples/color/colormaps_reference.html
    Returns :
    * color * rgba object, can be used in plot
    
    c = getcolor(10,0)
    ax.plot(x,y,color=c)  
    '''
    import matplotlib.pyplot as pl
    cmap = pl.cm.get_cmap(colormap)
    fact = float(itime)/ntimes
    rgba = cmap(fact)
    return rgba

def writetable(filename, M, rows_names_list,column_names_list,fmt="%5.3f\t"):
    '''
    Writes a 2d numpy array into a text file.
    It is a wrapper of np.savetxt

    Arguments:
    * filename * string
    * M        * a 2d numpy array
    * rows_names_list   * list of string, it will appear as first column
    * column_names_list * list of string, it will appear as first row
    
    
    Example:
    import numpy as np
    A=np.random.randn(3,2)
    writetable("tmp.txt", A, ["A1","A2","A3"], ['field1','field2'])
    writetable("tmp.txt", A, ["A1","A2","A3"], ['field1','field2'],fmt='%5.3f\t')
    writetable("tmp.txt", A, ["A1","A2","A3"], ['field1','field2'],fmt='%3.2f\t %d')
    '''

    headerstr="\t"
    nrows,ncols=M.shape
    dtype = [('colnames','U20')]
    if fmt.count("%") == 1:
        fmtlist = [fmt for ii in range(ncols)]
    else:
        fmtlist = fmt.split()
    typelist = []
    for ii in range(ncols):
        if 'f' in fmtlist[ii]: typelist.append(np.float32)
        if 's' in fmtlist[ii]: typelist.append('S15')
        if 'd' in fmtlist[ii]: typelist.append(int)
        
    for ii,col_name in enumerate(column_names_list):
        dtype.append((col_name,typelist[ii]))
        headerstr = headerstr + col_name + "\t "

    X = np.zeros((nrows,), dtype=dtype)
    for i in range(nrows):
        X['colnames'][i] = rows_names_list[i]
    for icol,colname  in enumerate(column_names_list):
        X[colname] = M[:,icol]
    if fmt.count("%") == 1:
        myformat = "%s\t" + fmt*ncols
    else:
        myformat = "%s\t" + fmt
    np.savetxt(filename, X, fmt=myformat, header=headerstr)

def data_for_linear_interp(array,value):
    '''
    Returns minimal data useful to perform linear interpolation
    two indexes and a weight.

    Arguments:
    * array * a sorted array, meaning x in a function y=f(x)
    * value * integer or float, meaning a point x(0)

    Returns :
    * before   * integer
    * after    * integer
    * t_interp * float


    value = array[before]*(1-t_interp) + array[after]*t_interp
    '''

    if value > array[-1] :
        l = len(array)
        return l-1, l-1, 0
    for i, array_elem in enumerate(array):
        if value < array_elem: break


    after = i
    before= i-1
    t_interp = float((value-array[before]))/(array[after]-array[before])

    if after==0 :
        before=0
        t_interp=0

    return before, after, t_interp

def Time_Interpolation(Instant_datetime, TimeList):
    '''
    Returns minimal data useful to perform linear interpolation in time,
    two indexes and a weight.

    Arguments:
    Instant_datetime : a datetime object
    Timelist         : a list of sorted datetime objects

    Returns :
    * before   * integer
    * after    * integer
    * t_interp * float

    '''
    instant_seconds = int(Instant_datetime.strftime("%s"))
    array___seconds = np.array([int(date.strftime("%s"))  for date in TimeList   ])

    before, after,t_interp = data_for_linear_interp(array___seconds,instant_seconds)
    return before, after,t_interp


def nan_compare(array_with_nan,operator, value):
    '''
    Compares array with scalar value, returning an array of logicals,
    without producing boresome warnings because of nans.

    Arguments:
    * array_with_nan * numpy array, it can contain nans
    * operator       * string, like ">", or ">="
    * value          * integer, or float

    nan_compare(array,">",value)  works like array > value
    It is the same of disable warnings, but it is safer.

    Returns:
    * out *  array of logicals, having the shape of array_with_nan.
    '''

    ii = np.isnan(array_with_nan)
    out = np.zeros(array_with_nan.shape, bool)
    if operator=='>'  : out[~ii] = array_with_nan[~ii] >  value
    if operator=='>=' : out[~ii] = array_with_nan[~ii] >= value
    if operator=='< ' : out[~ii] = array_with_nan[~ii] <  value
    if operator=='<=' : out[~ii] = array_with_nan[~ii] <= value
    return out


def nanmean_without_warnings(array):
    if len(array)==0: return np.nan
    if np.isnan(array).all(): return np.nan
    return np.nanmean(array)


def search_closest_sorted(a, b):
    """
    Returns an array `k` of the same shape of `b`. If `j` is the value of
    `k[i]`, this means that `a[j]` is the closest element (among all the
    elements of `a`) to `b[i]`.
    The vector `a` must be sorted (`a[i]` < `a[j]` for each `i` < `j`).

    In other words, this function behaves like `numpy.search_sorted` but,
    instead of returning always the element on the left (or on the right),
    it returns the closest one.

    :param a: 1D *sorted* numpy array
    :param b: numpy array or float
    """
    is_float = not hasattr(b, "__len__")
    b = np.asarray(b)

    a = np.asarray(a)
    if len(a.shape) != 1:
        raise ValueError('a must be a 1D array')

    if a.shape[0] == 0:
        raise ValueError('a can not be an empty array')

    on_left = b <= a[0]
    on_right = b > a[-1]
    inside_array = np.logical_not(np.logical_or(on_left, on_right))

    output = np.empty_like(b, dtype=np.int32)

    b_inside = b[inside_array]

    right_neighbour = np.searchsorted(a, b_inside, side="left")
    left_neighbour = right_neighbour - 1

    left_distances = np.abs(a[left_neighbour] - b_inside)
    right_distances = np.abs(a[right_neighbour] - b_inside)
    choose_left = left_distances < right_distances

    left_b_mask = np.zeros_like(b, dtype=bool)
    right_b_mask = np.zeros_like(b, dtype=bool)
    left_b_mask[inside_array] = choose_left
    right_b_mask[inside_array] = np.logical_not(choose_left)

    output[on_left] = 0
    output[on_right] = a.shape[0] - 1

    output[left_b_mask] = left_neighbour[choose_left]
    output[right_b_mask] = right_neighbour[np.logical_not(choose_left)]

    if is_float:
        return output[0]

    return output


if __name__ == '__main__':
    A=np.random.randn(3,2)
    writetable("tmp.txt", A, ["A1","A2","A3"], ['field1','field2'])
    writetable("tmp.txt", A, ["A1","A2","A3"], ['field1','field2'],fmt='%5.3f\t')

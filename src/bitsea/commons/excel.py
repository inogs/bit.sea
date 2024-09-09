import xlrd
import numpy as np

def read_sheet(filename, sheet_name):
    '''
    A, STRINGS = read_sheet(filename, sheet_name)
    Returns:
    * A       * numpy 2d array (nrows, ncols) of floats
    * STRINGS * numpy 2d array (nrows, ncols) of strings
    '''
    book = xlrd.open_workbook(filename)
    sheet = book.sheet_by_name(sheet_name)
    ncols = sheet.ncols
    nrows = sheet.nrows
    
    A       =  np.zeros((nrows,ncols), np.float32)
    STRINGS =  np.zeros((nrows,ncols), "S100")
    for j in range(ncols):
        COL=sheet.col(j)
        for i in range(nrows):
            t = COL[i]
            if t.ctype==2:
                A[i,j] = t.value
            else:
                STRINGS[i,j] = str(t.value)
               
    return A, STRINGS

if __name__ == "__main__":
    xls_file="/Users/gbolzon/Downloads/RIVERINE_INPUT_PERSEUS_BAU_D4.6_100.xls"
    A, scritte = read_sheet(xls_file, "monthly")
    B, scritte = read_sheet(xls_file, "DIP_KTperYR_NOBLS")
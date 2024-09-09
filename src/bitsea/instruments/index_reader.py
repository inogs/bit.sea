import numpy as np
from commons.utils import addsep
import os
def index_reader():
    GSS_DEFAULT_LOC = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/"
    ONLINE_REPO = addsep( os.getenv("ONLINE_REPO",GSS_DEFAULT_LOC))
    filename = ONLINE_REPO + "COPERNICUS/index_monthly.txt"
    mydtype= np.dtype([('catalog_id','S20'),
              ('file_name','S200'),
              ('geospatial_lat_min',np.float32),
              ('geospatial_lat_max',np.float32),
              ('geospatial_lon_min',np.float32),
              ('geospatial_lon_max',np.float32),
              ('time_coverage_start','S19'),
              ('time_coverage_end','S19'),
              ('provider','S30'),
              ('date_update','S30'),
              ('data_mode','S1'),
              ('parameters','S200')] )

    ALLVARLIST=['PHOS','SLCA','AMON','DOX1','DOX2','CPHL','NTRZ','NTRA','NTRI','PHPH']
    #NTRZ:long_name = "nitrate + nitrite" ;

    A=np.loadtxt(filename,dtype=mydtype, delimiter=",", skiprows=6 )
    return A


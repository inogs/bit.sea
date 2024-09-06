#! /usr/bin/python

import numpy as np

BIOOPTIMOD_type      = np.dtype([('PRES'       , np.float)   , ('VALUE'      , np.float)  ,  ('qc'    ,   'S100')   ,\
                                ('Latitude'    , np.float)   , ('Longitude'  , np.float)  ,  ('Julian_Day', np.float)   ,\
                                ('Name'        , 'S100')     , ('date'       , 'S100')])


QC_380_MEDSEA_MAY2019_BIOOPTIMOD_type      = np.dtype([('PRES'       , np.float)   , ('IRR_380'    , np.float)  ,  ('qc_380'    ,   'S100')   ,\
                                                       ('Latitude'   , np.float)   , ('Longitude'  , np.float)  ,  ('Julian_Day', np.float)   ,\
                                                       ('Name'       , 'S100')     , ('date'       , 'S100')])

QC_412_MEDSEA_MAY2019_BIOOPTIMOD_type      = np.dtype([('PRES'       , np.float)   , ('IRR_412'    , np.float)  ,  ('qc_412'    ,   'S100')   ,\
                                                       ('Latitude'   , np.float)   , ('Longitude'  , np.float)  ,  ('Julian_Day', np.float)   ,\
                                                       ('Name'       , 'S100')     , ('date'       , 'S100')])

QC_490_MEDSEA_MAY2019_BIOOPTIMOD_type      = np.dtype([('PRES'       , np.float)   , ('IRR_490'    , np.float)  ,  ('qc_490'    ,   'S100')   ,\
                                                       ('Latitude'   , np.float)   , ('Longitude'  , np.float)  ,  ('Julian_Day', np.float)   ,\
                                                       ('Name'       , 'S100')     , ('date'       , 'S100')])

QC_CHL_MEDSEA_MAY2019_BIOOPTIMOD_type      = np.dtype([('PRES'       , np.float)   , ('CHL'        , np.float)  ,  ('qc_chl'    ,   'S100')   ,\
                                                       ('Latitude'   , np.float)   , ('Longitude'  , np.float)  ,  ('Julian_Day', np.float)   ,\
                                                       ('Name'       , 'S100')     , ('date'       , 'S100')])

QC_PAR_MEDSEA_MAY2019_BIOOPTIMOD_type      = np.dtype([('PRES'       , np.float)   , ('PAR'        , np.float)  ,  ('qc_PAR'    ,   'S100')   ,\
                                                       ('Latitude'   , np.float)   , ('Longitude'  , np.float)  ,  ('Julian_Day', np.float)   ,\
                                                       ('Name'       , 'S100')     , ('date'       , 'S100')])

FloatIndex_type= np.dtype([
          ('file_name','S200'),
          ('lat',np.float32),
          ('lon',np.float32),
          ('time','S17'),
          ('parameters','S200')] )


plot_conf=np.dtype([('title'          ,'S100')  , ('var' ,'S100')        ,\
                    ('unit_conversion',np.float), ('unit_label','S100')  ,\
                    ('operator'       ,'S100')                           ,\
                    ('vmin'           ,np.float), ('vmax',np.float)      ,\
                    ('depth'          ,np.float)                         ,\
                    ('lonS'           ,np.float), ('lonE',np.float)      ,\
                    ('latS'           ,np.float), ('latE',np.float)      ,\
                    ('start'          ,'S100')  , ('end'   ,'S100')      ,\
                    ('indataAve'      ,'S100')  , ('aggregation' ,'S100'),\
                    ('InputDir'       ,'S100')  , ('OutputDir'   ,'S100'),\
                    ('FileName'       ,'S100')]) 


Ed380_type=    [ ('profile'       ,'S100') ,    ('depth',     np.float32),   ('IRR_380',     np.float32)     ,\
                 ('latitude',  np.float32) ,    ('longitude', np.float32),   ('year',                np.int32)     ,\
                 ('month',        np.int32),    ('day',         np.int32),   ('basin',                 'S100') ]

Ed412_type=    [ ('profile'       ,'S100') ,    ('depth',     np.float32),   ('IRR_412',     np.float32)     ,\
                 ('latitude',  np.float32) ,    ('longitude', np.float32),   ('year',                np.int32)     ,\
                 ('month',        np.int32),    ('day',         np.int32),   ('basin',                 'S100') ]

Ed490_type=    [ ('profile'       ,'S100') ,    ('depth',     np.float32),   ('IRR_490',     np.float32)     ,\
                 ('latitude',  np.float32) ,    ('longitude', np.float32),   ('year',                np.int32)     ,\
                 ('month',        np.int32),    ('day',         np.int32),   ('basin',                 'S100') ]

PAR_type=    [ ('profile'       ,'S100') ,    ('depth',     np.float32),   ('PAR',     np.float32)     ,\
                 ('latitude',  np.float32) ,    ('longitude', np.float32),   ('year',                np.int32)     ,\
                 ('month',        np.int32),    ('day',         np.int32),   ('basin',                 'S100') ]
#! /usr/bin/python

import numpy as np

CHLOROPHYLL_PROFILE_QC_MEDSEA_type  = np.dtype([('name_profile' , 'S100')    , ('PRES'           , np.float)                            ,\
                                             ('JULD'            , np.float)  , ('CHLA_CALIBRATED', np.float)                            ,\
                                             ('Date'            , 'S100')    , ('cycle'          , np.int)    , ('profile', np.int)     ,\
                                             ('lat'             , np.float)  , ('long'           , np.float)                            ,\
                                             ('bon'             , np.int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , np.int)    , ('month'          , np.int)    , ('day'    , np.int) ])

Ed380_PROFILE_TYPE_1_MEDSEA_type    = np.dtype([('name_profile' , 'S100')    , ('PRES'           , np.float)                            ,\
                                             ('IRR_380'         , np.float)  , ('qc_380'         , 'S100')                              ,\
                                             ('Date'            , 'S100')    , ('cycle'          , np.int)    , ('profile', np.int)     ,\
                                             ('lat'             , np.float)  , ('long'           , np.float)                            ,\
                                             ('bon'             , np.int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , np.int)    , ('month'          , np.int)    , ('day'    , np.int) ])

 
Ed412_PROFILE_TYPE_1_MEDSEA_type    = np.dtype([('name_profile' , 'S100')    , ('PRES'           , np.float)                            ,\
                                             ('IRR_412'         , np.float)  , ('qc_412'         , 'S100')                              ,\
                                             ('Date'            , 'S100')    , ('cycle'          , np.int)    , ('profile', np.int)     ,\
                                             ('lat'             , np.float)  , ('long'           , np.float)                            ,\
                                             ('bon'             , np.int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , np.int)    , ('month'          , np.int)    , ('day'    , np.int) ])

Ed490_PROFILE_TYPE_1_MEDSEA_type    = np.dtype([('name_profile' , 'S100')    , ('PRES'           , np.float)                            ,\
                                             ('IRR_490'         , np.float)  , ('qc_490'         , 'S100')                              ,\
                                             ('Date'            , 'S100')    , ('cycle'          , np.int)    , ('profile', np.int)     ,\
                                             ('lat'             , np.float)  , ('long'           , np.float)                            ,\
                                             ('bon'             , np.int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , np.int)    , ('month'          , np.int)    , ('day'    , np.int) ])

PAR_PROFILE_TYPE_1_MEDSEA_type      = np.dtype([('name_profile' , 'S100')    , ('PRES'           , np.float)                            ,\
                                             ('IRR_PAR'         , np.float)  , ('qc_PAR'         , np.float)                            ,\
                                             ('Date'            , 'S100')    , ('cycle'          , np.int)    , ('profile', np.int)     ,\
                                             ('lat'             , np.float)  , ('long'           , np.float)                            ,\
                                             ('bon'             , np.int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , np.int)    , ('month'          , np.int)    , ('day'    , np.int) ])

SALINITY_PROFILE_MEDSEA_type        = np.dtype([('name_profile' , 'S100')    , ('PRES'           , np.float)                            ,\
                                             ('JULD'            , np.float)  , ('SAL'            , np.float)                            ,\
                                             ('Date'            , 'S100')    , ('cycle'          , np.int)    , ('profile', np.int)     ,\
                                             ('lat'             , np.float)  , ('long'           , np.float)                            ,\
                                             ('bon'             , np.int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , np.int)    , ('month'          , np.int)    , ('day'    , np.int) ])

TEMPERATURE_PROFILE_MEDSEA_type     = np.dtype([('name_profile' , 'S100')    , ('PRES'           , np.float)                            ,\
                                             ('JULD'            , np.float)  , ('TEMP'           , np.float)                            ,\
                                             ('Date'            , 'S100')    , ('cycle'          , np.int)    , ('profile', np.int)     ,\
                                             ('lat'             , np.float)  , ('long'           , np.float)                            ,\
                                             ('bon'             , np.int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , np.int)    , ('month'          , np.int)    , ('day'    , np.int) ])
List_QC_profiles_MedSea_type        = np.dtype([('name_profile' , 'S100')    , ('Zeu'            , np.float)                            ,\
                                             ('Zpd'             , np.float)  , ('Ed_inw_PAR'     , np.float)                            ,\
                                             ('Ed_inw_380'      , np.float)  , ('Ed_inw_412'     , np.float)                            ,\
                                             ('Ed_inw_490'      , np.float)  , ('Date'           , 'S100')                              ,\
                                             ('cycle'           , np.int)    , ('profile'        , np.int)                              ,\
                                             ('lat'             , np.float)  , ('long'           , np.float)                            ,\
                                             ('bon'             , np.int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , np.int)    , ('month'          , np.int)    , ('day'    , np.int) ])


FloatIndex_type= np.dtype([
          ('file_name','S200'),
          ('lat',np.float32),
          ('lon',np.float32),
          ('time','S17'),
          ('parameters','S200')] )

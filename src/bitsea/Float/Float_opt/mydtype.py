#! /usr/bin/python

import numpy as np

CHLOROPHYLL_PROFILE_QC_MEDSEA_type  = np.dtype([('name_profile' , 'S100')    , ('PRES'           , float)                            ,\
                                             ('JULD'            , float)  , ('CHLA_CALIBRATED', float)                            ,\
                                             ('Date'            , 'S100')    , ('cycle'          , int)    , ('profile', int)     ,\
                                             ('lat'             , float)  , ('long'           , float)                            ,\
                                             ('bon'             , int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , int)    , ('month'          , int)    , ('day'    , int) ])

Ed380_PROFILE_TYPE_1_MEDSEA_type    = np.dtype([('name_profile' , 'S100')    , ('PRES'           , float)                            ,\
                                             ('IRR_380'         , float)  , ('qc_380'         , 'S100')                              ,\
                                             ('Date'            , 'S100')    , ('cycle'          , int)    , ('profile', int)     ,\
                                             ('lat'             , float)  , ('long'           , float)                            ,\
                                             ('bon'             , int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , int)    , ('month'          , int)    , ('day'    , int) ])

 
Ed412_PROFILE_TYPE_1_MEDSEA_type    = np.dtype([('name_profile' , 'S100')    , ('PRES'           , float)                            ,\
                                             ('IRR_412'         , float)  , ('qc_412'         , 'S100')                              ,\
                                             ('Date'            , 'S100')    , ('cycle'          , int)    , ('profile', int)     ,\
                                             ('lat'             , float)  , ('long'           , float)                            ,\
                                             ('bon'             , int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , int)    , ('month'          , int)    , ('day'    , int) ])

Ed490_PROFILE_TYPE_1_MEDSEA_type    = np.dtype([('name_profile' , 'S100')    , ('PRES'           , float)                            ,\
                                             ('IRR_490'         , float)  , ('qc_490'         , 'S100')                              ,\
                                             ('Date'            , 'S100')    , ('cycle'          , int)    , ('profile', int)     ,\
                                             ('lat'             , float)  , ('long'           , float)                            ,\
                                             ('bon'             , int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , int)    , ('month'          , int)    , ('day'    , int) ])

PAR_PROFILE_TYPE_1_MEDSEA_type      = np.dtype([('name_profile' , 'S100')    , ('PRES'           , float)                            ,\
                                             ('IRR_PAR'         , float)  , ('qc_PAR'         , float)                            ,\
                                             ('Date'            , 'S100')    , ('cycle'          , int)    , ('profile', int)     ,\
                                             ('lat'             , float)  , ('long'           , float)                            ,\
                                             ('bon'             , int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , int)    , ('month'          , int)    , ('day'    , int) ])

SALINITY_PROFILE_MEDSEA_type        = np.dtype([('name_profile' , 'S100')    , ('PRES'           , float)                            ,\
                                             ('JULD'            , float)  , ('SAL'            , float)                            ,\
                                             ('Date'            , 'S100')    , ('cycle'          , int)    , ('profile', int)     ,\
                                             ('lat'             , float)  , ('long'           , float)                            ,\
                                             ('bon'             , int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , int)    , ('month'          , int)    , ('day'    , int) ])

TEMPERATURE_PROFILE_MEDSEA_type     = np.dtype([('name_profile' , 'S100')    , ('PRES'           , float)                            ,\
                                             ('JULD'            , float)  , ('TEMP'           , float)                            ,\
                                             ('Date'            , 'S100')    , ('cycle'          , int)    , ('profile', int)     ,\
                                             ('lat'             , float)  , ('long'           , float)                            ,\
                                             ('bon'             , int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , int)    , ('month'          , int)    , ('day'    , int) ])
List_QC_profiles_MedSea_type        = np.dtype([('name_profile' , 'S100')    , ('Zeu'            , float)                            ,\
                                             ('Zpd'             , float)  , ('Ed_inw_PAR'     , float)                            ,\
                                             ('Ed_inw_380'      , float)  , ('Ed_inw_412'     , float)                            ,\
                                             ('Ed_inw_490'      , float)  , ('Date'           , 'S100')                              ,\
                                             ('cycle'           , int)    , ('profile'        , int)                              ,\
                                             ('lat'             , float)  , ('long'           , float)                            ,\
                                             ('bon'             , int)    , ('float'          , 'S100')                              ,\
                                             ('project'         , 'S100')    , ('type'           , 'S100')                              ,\
                                             ('basin'           , 'S100')    , ('region'         , 'S100')                              ,\
                                             ('year'            , int)    , ('month'          , int)    , ('day'    , int) ])


FloatIndex_type= np.dtype([
          ('file_name','S200'),
          ('lat',np.float32),
          ('lon',np.float32),
          ('time','S17'),
          ('parameters','S200')] )

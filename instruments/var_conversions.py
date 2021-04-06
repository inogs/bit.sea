FLOATVARS={'O2o':'DOXY', \
           'N3n':'NITRATE',  \
           'P_l':'CHLA', \
           'P_i':'CHLA', \
           'Chla':'CHLA', \
           'vosaline':'PSAL', \
           'votemper':'TEMP',
           'EIR':'DOWNWELLING_PAR',
           'POC':'BBP700',
           'P_c':'BBP700',
           'BBP532':'BBP532',
           'pH'    : 'PH_IN_SITU_TOTAL'   }
 
#DOXY:units = "micromole/kg" ;
#CHLA:units = "mg/m3" ;
#NITRATE:units = "micromole/kg" ;

LOVFLOATVARS={'O2o':'DOXY', \
              'N3n':'SR_NO3',  \
              'P_l':'CHLA', \
              'Chla':'CHLA', \
              'P_i':'CHLA', \
              'vosaline':'PSAL', \
              'votemper':'TEMP',
              'EIR':'PAR',
              'BBP532':'BBP532',
              'POC':'BBP700',
              'P_c':'BBP700',
               'pH'    : 'PH_IN_SITU_TOTAL'  }



MOORINGVARS={'O2o':'DOX1', \
             'N3n':'NOTFOUND',  \
             'P_i':'CPHL' , \
             'P_l':'CPHL'}
#DOX1:units = "ml/l" ;
#CPHL:units = "milligram/m3" ;


VESSELVARS={'O2o':'DOX1', \
             'N3n':'NTRA', \
             'P_i':'CPHL', \
             'P_l':'CPHL', \
             'N1p':'PHOS', \
             'N5s':'SLCA', \
             'pH_':'PHPH'}

CARBONVARS={'O2o':'oxygen', \
            'N3n':'nitrate',  \
            'N1p':'phosphate', \
            'N5s':'silicate', \
            'ALK':'ALK', \
            'DIC':'DICric', \
            'pH' :'PHt_{T-Press-ins}', \
           'pCO2':'pCO2'}

NUTRVARS={'O2o':'oxygen', \
             'N3n':'nitrate',  \
             'N1p':'phosphate', \
             'N4n':'ammonium', \
             'P_l': 'chlorophyll',\
             'N5s':'silicate'}


ISPRAVARS = {'N1p': 'orthophosphates', \
             'N3n': 'nitrate', \
             'P_l': 'chlorophyll a'}

MASSIMILIVARS = {'N1p': 'phosphate', \
             'N3n': 'nitrate', \
             'P_l': 'CHL'}

FLOAT_OPT_VARS = {"P_l": "chl",
                  "Ed_380": "Ed_380",
                  "Ed_412": "Ed_412",
                  "Ed_490": "Ed_490",
                  }

FLOAT_OPT_VARS_2019 = {"P_l": "CHL",
                  "Ed_380": "IRR_380",
                  "Ed_412": "IRR_412",
                  "Ed_490": "IRR_490",
                  "Ed_par": "PAR",
                  }

SUPERFLOAT_VARS = {"P_l"    : "CHLA",
                  "Ed_380"  : "IRR_380",
                  "Ed_412"  : "IRR_412",
                  "Ed_490"  : "IRR_490",
                  "Ed_par"  : "PAR",
                  "O2o"     : "DOXY",
                  "N3n"     : "NITRATE",
                  "CDOM"    : "CDOM",
                  "vosaline": "SALI",
                  "votemper": "TEMP",
                  "POC"     : "BBP700",
                  }

SOCAT_VARS ={"votemper": "temp",
             "vosaline": "sal",
             "pCO2"    : "fCO2"
             }



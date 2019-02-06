FLOATVARS={'O2o':'DOXY', \
           'N3n':'NITRATE',  \
           'P_l':'CHLA', \
           'P_i':'CHLA', \
           'Chla':'CHLA', \
           'vosaline':'PSAL', \
           'votemper':'TEMP' }   
#DOXY:units = "micromole/kg" ;
#CHLA:units = "mg/m3" ;
#NITRATE:units = "micromole/kg" ;

LOVFLOATVARS={'O2o':'DOXY', \
              'N3n':'SR_NO3',  \
              'P_l':'CHLA', \
              'Chla':'CHLA', \
              'P_i':'CHLA', \
              'vosaline':'PSAL', \
              'votemper':'TEMP' }



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
            'N5s':'silicate'}

NUTRVARS={'O2o':'oxygen', \
             'N3n':'nitrate',  \
             'N1p':'phosphate', \
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

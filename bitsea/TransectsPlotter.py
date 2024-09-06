#!/usr/bin/env python
# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

from __future__ import print_function

import sys
import os
import re
from glob import glob
from warnings import warn

try:
    from transects.transect import Transect
    from transects.plot import transectplot
    from commons.mask import Mask
    from commons.dataextractor import NotFoundError
except ImportError:
    print("You should run this script from the bit.sea root directory.", file=sys.stderr)
    sys.exit(2)

#Input directory
INPUTDIR = 'layer_integral'

#Input file glob
INPUTGLOB = 'ave*nc'

#Mask file
MESHMASK = 'layer_integral/meshmask.nc'

#Plot list XML file
PLOTLISTFILE = 'postproc/Plotlist.xml'

#Output directory
OUTPUTDIR = '/tmp/output'

if __name__ == "__main__":
    #Build transect list from Plotlist.xml
    try:
        transect_list = Transect.get_list_from_XML_file(PLOTLISTFILE)
    except Exception as e:
        msg = "%s:\n%s" % (PLOTLISTFILE, e)
        raise Exception(msg)

    #Open mask file
    try:
        mask = Mask(MESHMASK)
    except Exception as e:
        msg = "Unable to open mask file %s:\n%s" % (meshmask, e)
        raise Exception(msg)

    #Build input files list
    inputfiles = glob(os.path.join(INPUTDIR,INPUTGLOB))

    #Check if output directory exists
    if os.path.exists(OUTPUTDIR):
        #Verify that it is a directory
        if not os.path.isdir(OUTPUTDIR):
            raise Exception("'%s' is not a directory." % str(OUTPUTDIR))
    else:
        #Try to create it
        os.mkdir(OUTPUTDIR)

    #For each input file
    for nc_file in inputfiles:
        #Get the date
        m = re.search('.*([0-9]{4})([0-9]{2})([0-9]{2}).*',nc_file)
        date = "%s-%s-%s" % (m.group(1), m.group(2), m.group(3))
        shortdate = "%s%s%s" % (m.group(1), m.group(2), m.group(3))
        #For each Transect
        for t in transect_list:
            #Retrieve data from file
            try:
                t.fill_data_cache_from_file(mask, nc_file)
            except NotFoundError as e:
                warn(str(e))
                continue
            sd = None
            #For each segment
            for s in t.segmentlist:
                #Plot segment data
                sd, fig, ax = transectplot(t, s, date, segmentdata=sd, dpi=18)
                #Save figure
                outfilename = "%s/ave.%s.%s.%s.png" % (OUTPUTDIR, shortdate, t.varname, s.name)
                fig.savefig(outfilename)

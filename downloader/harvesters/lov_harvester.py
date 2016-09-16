import urllib2
import os
from vlfr import username, password

web_site="http://www.oao.obs-vlfr.fr/BD_FLOAT/NETCDF/"
OUTPUTDIR = "/Users/gbolzon/Downloads/FLOATS/"

manager = urllib2.HTTPPasswordMgrWithDefaultRealm()
manager.add_password(None, web_site, username, password)

auth = urllib2.HTTPBasicAuthHandler(manager)
opener = urllib2.build_opener(auth)
urllib2.install_opener(opener)




float_name="lovbio085d"
REMOTEDIR=web_site + float_name


urlfilelist = REMOTEDIR + "/liste_all"
response = urllib2.urlopen(urlfilelist)
remotefilelist = response.read().rsplit("\n")[:-1]

filelist=[os.path.basename(f) for f in remotefilelist]


for filename in filelist:
    localfile = OUTPUTDIR + filename
    if os.path.exists(localfile):
        continue
    print "downloading " , filename
    url = REMOTEDIR + "/" + filename    
    response = urllib2.urlopen(url)
    F = open(localfile,'wb')
    F.write(response.read())
    F.close()
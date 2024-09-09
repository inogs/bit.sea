# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
from __future__ import print_function

import sys
import os,code
import urllib2
from os.path import join, exists, realpath, dirname
from os import remove, rename
from xml.dom.minidom import parseString

import xml.etree.ElementTree as xml_tree

from harvester_interface import HarvesterInterface
from utilities.files_and_dirs import ensure_dir
from utilities.date_and_time import now_as_string
from vlfr import username, password
import numpy as np

http_url = "http://www.oao.obs-vlfr.fr/BD_FLOAT/NETCDF/"
relative_path = "FLOAT_LOVBIO"

manager = urllib2.HTTPPasswordMgrWithDefaultRealm()        
manager.add_password(None, http_url, username, password)

auth = urllib2.HTTPBasicAuthHandler(manager)
opener = urllib2.build_opener(auth)
urllib2.install_opener(opener)

wmo_default = realpath(dirname(realpath(__file__)) + "/../harvesters_info/wmo.txt")
wmo_file=os.getenv("WMO_FILE")
if wmo_file is None:
    print("Env WMO_FILE is not defined. Taking the hardcoded wmo.txt")
    wmo_file=wmo_default
else:
    print("Taking regularly " + wmo_file)


xml_path = realpath(dirname(realpath(__file__)) + '/../harvesters_xml')
print(wmo_file)

class ErrorMessage(object):
    """A class to handle error messages"""
    def __init__(self):
        # self.message contains the error
        self.message = None
def download_file(url, f, path, log, perms=None,
                  skip_if_exists=True, skip_is_strange = False):
    """
    Download a file from an ftp remote archive.
    
    Args:
        - *connection*: An ftp connection, created with the ftplib library
        - *f*: The name of the file as a string
        - *path*: The directory where the file will be saved
        - *log*: A Log object to report the operations to the user
        - *perms* (optional): A dictionary of permissions like the one the
          list_file function returns. If this dictionary is not none, the
          file will be downloaded only if the readable flag for others is
          present.
        - *skip_if_exists*: A boolean value that says what to do if the file
          is already present. Default is True.
        - *skip_is_strange*: This means that the fact that the file is already
          present shall be logged with high severity. If skip_if_exists is
          False it has no effect. Default is False.
    
    Returns:
        *Success*: A boolean value that reports if the downloade
        was executed or if it was skipped for some reasons
    
    Raises:
        *DownloadError* if something went wrong during the download
    """

    outfile = join(path, f)
    if skip_if_exists:
        # Check if the file already exists
        if exists(outfile):
            if skip_is_strange:
                log.info('Skipping file ' + f + ' because'
                         ' it was already downloaded!')
            else:
                log.debug('Skipping file ' + f + ' because '
                          'it was already downloaded')
            return False

    # Check if the file is readable
    if perms is not None:
        if 'read' not in perms[f]['others']:
            log.info('The file ' + f + ' is not '
                     'readable. Check the '
                     'permissions')
            return False
 
    # Download the file
    log.debug('Downloading ' + f + '...')
    
    # I will first save the partial file with a temp name. When the
    # transfers end, I will rename it

    temp_name = join(path, 'incomplete_download.tmp')
    response = urllib2.urlopen(url)
    F = open(temp_name,'wb')
    F.write(response.read())
    F.close()

    # Move the temp file in the new downloaded file
    rename(temp_name, outfile)
    
    log.debug('File ' + f + ' downloaded!')
    return True



class LovBioFloatsHarvester(HarvesterInterface):
    """
    This is the harvester in charge of download all the files whose name
    start with "MR" from the ftp server ftp.ifremer.fr. This harvester
    need a file (called wmo_file) that reports the status of the biological
    floats. The position of that file is set in the global variable "wmo_file".
    Moreover, the harvester will save the information of the previous file
    in a xml file whose path is set in the global variable "xml_path".
    """
    def is_a_lov_float(self,wmo,lovname):
        if wmo == lovname : return False
        if lovname == '' : return False
        return True

    def wmo_file_reader(self):
        FLOAT_dtype=[('id',np.int),('wmo','S20'),('id_type',np.int),('type','S20'),('nome_fs','S20'),('status','S1')]
        TABLE=np.loadtxt(wmo_file,dtype=FLOAT_dtype,skiprows=2,delimiter=" | ")
        return TABLE

    def harvest(self, db_path, log):
        """
        For every float in the file wmo, check the status in the wmo_file and
        in the xml one. If in at least one file the float is reported as active,
        then check the last file downloaded for that wmo and download every file
        on the server that is more recent than the one already downloaded. Then
        update the xml file with the status reported in the wmo file.

        Args:
            - *db_path*: the path of the directory set in the download program.
            - *log*: a logger object from the class Log to print informations on
              the standard output
               
        Returns:
            - *downloaded*: a list of all the downloaded filenames.
        """
        # In the following list I will store the name of the
        # files that will be downloaded or updated
        print("HARVEST")
        downloaded = []

        # Read the wmo file

        A = self.wmo_file_reader()
        lines_active_floats=np.where(A['status']=='A')[0]
        lines_dead__floats =np.where(A['status']=='D')[0]


        # Now we need the xml file that keeps what we did on the
        # last updates
        xml_file = join(xml_path, self.__class__.__name__ + '.xml')
        try:
            tree = xml_tree.parse(xml_file)
        except:
            log.info('XML file not found or not readable. '
                     'This script will update every file '
                     'from the remote archive. This is '
                     'almost the same than run in reset '
                     'mode, but the files that exist will '
                     'not be downloaded again. Moreover, '
                     'the XML file will be rewritten.')
            return self.rebuild(db_path, log, skip_if_present=True)

        root = tree.getroot()

        # Check if the directory for this harvester is present
        # in the database
        path = join(db_path,relative_path)
        ensure_dir(path, log, expected=True)

        # Open the connection with the remote archive

        
        # Enter in the directory tree
        

        # Download data for every active float
        for l in lines_active_floats:
            # Update the xml with the current status of the float
            f = A[l]['wmo']
            floatname = A[l]['nome_fs'].replace(' ','')
            if not self.is_a_lov_float(f, floatname): continue

            wmo_in_xml = 'wmo_' + str(f)
            f_in_xml = root.findall(wmo_in_xml)
            if len(f_in_xml) == 0:
                f_node = xml_tree.SubElement(root, wmo_in_xml)
            else:
                f_node = [fn for fn in f_in_xml if fn.tag==wmo_in_xml][0]
            f_node.set('status', 'A')

            try:
                urlfilelist = http_url + floatname +  "/liste_all"
                print(urlfilelist)
                response = urllib2.urlopen(urlfilelist)
            except:
                log.info('No directory associated with file ' + str(f) +
                         '. This file will be skipped!')
                continue

            remotepathlist = response.read().rsplit("\n")[:-1]
            filelist=[os.path.basename(fn) for fn in remotepathlist]
            # Now I look for the profiles dir. This is the folder
            # where all the data are stored


            if len(filelist) > 0:
                download_for_f = []
                # Copy all file in a local dir with the same name
                # skipping the one that we already have
                float_local_dir = join(path, f)
                ensure_dir(float_local_dir, log, expected = False)
                for ff in filelist:
                    url = http_url + floatname + "/" + ff
                    d = download_file(url, ff, float_local_dir,
                                      log, None, True)
                    # If the file was downloaded without any problem,
                    # add it to the list of downloaded files
                    if d:
                        downloaded.append(ff)
                        download_for_f.append(ff)
                if len(download_for_f) == 0:
                    log.info('No updates found for float ' + str(f))                    
            else:
                log.info('No updates found for float ' + str(f))


        print ("DIED FLOATS")
        for l in lines_dead__floats:
            f = A[l]['wmo']
            floatname = A[l]['nome_fs'].replace(' ','')
            if not self.is_a_lov_float(f, floatname): continue

            to_be_updated = False
            # Update the xml with the current status of the float
            # Check if it must be updated
            f_in_xml = root.findall('wmo_' + str(f))
            if len(f_in_xml) == 0:
                # If this float is new, then add it to the archive
                # and it will be updated
                to_be_updated = True
                f_node = xml_tree.SubElement(root, 'wmo_' + str(f))
            else:
                f_node = [fn for fn in f_in_xml if fn.tag=='wmo_'+str(f)][0]
                # If I already know this float, but the last time it
                # was not dead, update it
                if f_node.get('status') != 'D':
                    to_be_updated = True
            f_node.set('status', 'D')
            
            if not to_be_updated:
                log.debug("Wmo " + str(f) + " is dead and will not be updated")
            else:
                log.debug("Wmo " + str(f) + " now is dead but was active on "
                          "the last run and will be updated anyway")


            if to_be_updated:
                try:
                    urlfilelist = http_url + floatname +  "/liste_all"
                    print(urlfilelist)
                    response = urllib2.urlopen(urlfilelist)
                except:
                    log.info('No directory associated with file ' + str(f) +
                             '. This file will be skipped!')
                    continue

                remotepathlist = response.read().rsplit("\n")[:-1]
                filelist=[os.path.basename(fn) for fn in remotepathlist]
                # Now I look for the profiles dir. This is the folder
                # where all the data are stored
                if len(filelist) > 0:
                    download_for_f = []
                    # Copy all file in a local dir with the same name
                    # skipping the one that we already have
                    float_local_dir = join(path, f)
                    ensure_dir(float_local_dir, log, expected = False)
                    for ff in filelist:
                        url = http_url + floatname + "/" + ff
                        d = download_file(url, ff, float_local_dir,
                                          log, None, True)
                        # If the file was downloaded without any problem,
                        # add it to the list of downloaded files
                        if d:
                            downloaded.append(ff)
                    if len(download_for_f) == 0:
                        log.info('No updates found for float ' + str(f))                    
                else:
                    log.info('No updates found for float ' + str(f))

        # Save the XML file
        root.set('Updated', now_as_string())
        xml_as_string = xml_tree.tostring(root)
        xml_rebuild = parseString(xml_as_string)
        pretty_xml = xml_rebuild.toprettyxml(indent='  ')
        pretty_xml_lines = pretty_xml.split('\n')
        pretty_xml = "\n".join([l for l in pretty_xml_lines if l.strip()])
        with open(xml_file, 'w') as xml_f:
            xml_f.write(pretty_xml)

        # Return the list of downloaded files
        return downloaded


    def rebuild(self, db_path, log, skip_if_present=False):
        """
        For every float in the file wmo, download every data file related to
        that float that starts with 'MR'. Then create a xml file with the
        data read from the wmo file.

        Args:
            - *db_path*: the path of the directory set in the download program.
            - *log*: a logger object from the class Log to print informations on
              the standard output
            - *skip_if_present*: A boolean value that set if not download again 
              the files that are saved on the local directory. By defalut is False
               
        Returns:
            - *downloaded*: a list of all the downloaded filenames.
        """
        # In the following list I will store the name of the
        # files that will be downloaded or updated
        print("REBUILD")
        downloaded = []

        # Read the wmo file line by line (exclude the first one because
        # it does not contain data)
        A=self.wmo_file_reader()


        # Delete, if present, the XML files with all the floats
        xml_file = join(xml_path, self.__class__.__name__ + '.xml')
        if exists(xml_file):
            remove(xml_file)
        # and create a new one (in memory)
        root = xml_tree.Element("BioFloats")
        root.set('Updated', now_as_string())
        tree = xml_tree.ElementTree(root)

        # Check if the directory for this harvester is present
        # in the database
        path = join(db_path,relative_path)
        ensure_dir(path, log, expected=False)

        # Open the connection with the remote archive

        # Download data for every active float
        for l in range(len(A)):
            f = A[l]['wmo']
            floatname = A[l]['nome_fs'].replace(' ','')
            if not self.is_a_lov_float(f, floatname): continue


            # Update the xml with the current status of the float
            f_in_xml = root.findall('wmo_' + str(f))
            if len(f_in_xml) == 0:
                f_node = xml_tree.SubElement(root, 'wmo_' + str(f))

            else:
                f_node = [fn for fn in f_in_xml if fn.tag=='wmo_'+str(f)][0]
            f_node.set('status', A[l]['status'])

            try:
                urlfilelist = http_url + floatname +  "/liste_all"
                print(urlfilelist)
                response = urllib2.urlopen(urlfilelist)
            except:
                log.info('Cannot download file ' + urlfilelist +
                         '. This file will be skipped!')
                continue

            remotepathlist = response.read().rsplit("\n")[:-1]
            filelist=[os.path.basename(fn) for fn in remotepathlist]
            # Now I look for the profiles dir. This is the folder
            # where all the data are stored
            if len(filelist) > 0:
                download_for_f = []
                # Copy all file in a local dir with the same name
                # skipping the one that we already have
                float_local_dir = join(path, f)

                print(float_local_dir)
                ensure_dir(float_local_dir, log, expected = False)
                for ff in filelist:
                    url = http_url + floatname + "/" + ff
                    d = download_file(url, ff, float_local_dir,
                                      log, None, True)
                    # If the file was downloaded without any problem,
                    # add it to the list of downloaded files
                    if d:
                        downloaded.append(ff)
                        download_for_f.append(ff)
                if len(download_for_f) == 0:
                    log.info('No updates found for float ' + str(f))                    
            else:
                log.info('No updates found for float ' + str(f))



        
        # Save the XML file
        xml_as_string = xml_tree.tostring(root)
        xml_rebuild = parseString(xml_as_string)
        pretty_xml = xml_rebuild.toprettyxml(indent='  ')
        pretty_xml_lines = pretty_xml.split('\n')
        pretty_xml = "\n".join([l for l in pretty_xml_lines if l.strip()])

        ensure_dir(xml_path, log, expected=False)
        with open(xml_file, 'w') as xml_f:
            xml_f.write(pretty_xml)

        # Return the list of downloaded files
        return downloaded

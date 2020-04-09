# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
from __future__ import print_function

import sys,os

from ftplib import FTP
from os.path import join, exists, realpath, dirname
from os import listdir, remove
from xml.dom.minidom import parseString

import xml.etree.ElementTree as xml_tree

from harvester_interface import HarvesterInterface
from utilities.ftp_utilities import list_files, download_file
from utilities.files_and_dirs import ensure_dir
from utilities.date_and_time import now_as_string
import numpy as np
ftp_url = 'ftp.ifremer.fr'

relative_path = "FLOAT_BIO"
wmo_default = realpath(dirname(realpath(__file__)) + "/../harvesters_info/wmo.txt")
wmo_file=os.getenv("WMO_FILE")
if wmo_file is None:
    print("Env WMO_FILE is not defined. Taking the hardcoded wmo.txt")
    wmo_file=wmo_default
else:
    print("Taking regularly " + wmo_file)


xml_path = realpath(dirname(realpath(__file__)) + '/../harvesters_xml')
print(wmo_file)


class BioFloatsHarvester(HarvesterInterface):
    """
    This is the harvester in charge of download all the files whose name
    start with "SR" or "SD" from the ftp server ftp.ifremer.fr. This harvester
    need a file (called wmo_file) that reports the status of the biological
    floats. The position of that file is set in the global variable "wmo_file".
    Moreover, the harvester will save the information of the previous file
    in a xml file whose path is set in the global variable "xml_path".
    """
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
        downloaded = []
        A = self.wmo_file_reader()
        #lines_active_floats=np.where(A['status']=='A')[0]
        #lines_dead__floats =np.where(A['status']=='D')[0]

        # Read the wmo file line by line (exclude the first one because
        # it does not contain data)
        #wmo_list = A['wmo']
        nFloats=A.size

        # Put its content in a dictionary
        wmo_status = dict()
        for l in range(nFloats):
            wmo_status[A[l]['wmo']] = A[l]['status']

        active_floats = [f for f in wmo_status if wmo_status[f]=='A']
        dead_floats = [f for f in wmo_status if wmo_status[f]=='D']

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
        connection = FTP(ftp_url)
        connection.login()

        # Enter in the directory tree
        connection.cwd('ifremer/argo/dac/coriolis')

        # Download data for every active float
        for f in active_floats:
            # Update the xml with the current status of the float
            f_in_xml = root.findall('wmo_' + str(f))
            if len(f_in_xml) == 0:
                f_node = xml_tree.SubElement(root, 'wmo_' + str(f))
                f_node.set('status', 'A')
            else:
                f_node = [fn for fn in f_in_xml if fn.tag=='wmo_'+str(f)][0]
                f_node.set('status', 'A')

            try:
                connection.cwd(f)
            except:
                log.info('No directory associated with wmo ' + str(f) +
                         '. This wmo will be skipped!')
                continue

            _, float_dirs, _ = list_files(connection)
            print(float_dirs)
            # Now I look for the profiles dir. This is the folder
            # where all the data are stored
            if 'profiles' in float_dirs:
                connection.cwd('profiles')
                float_files, _, perms = list_files(connection)
                min_len = len(f) + 2
                to_be_downloaded = [ff for ff in float_files
                                    if len(ff)>min_len
                                    and ff[:min_len] in ['SR'+f,'SD'+f ]]
                if len(to_be_downloaded) > 0:
                    download_for_f = []
                    # Copy all file in a local dir with the same name
                    # skipping the one that we already have
                    float_local_dir = join(path, f)
                    ensure_dir(float_local_dir, log, expected = False)
                    for ff in to_be_downloaded:
                        d = download_file(connection, ff, float_local_dir,
                                          log, perms, True)
                        # If the file was downloaded without any problem,
                        # add it to the list of downloaded files
                        if d:
                            downloaded.append(ff)
                            download_for_f.append(ff)
                    if len(download_for_f) == 0:
                        log.info('No updates found for float ' + str(f))                    
                else:
                    log.info('No updates found for float ' + str(f))
                connection.cwd('..')
            else:
                log.info('Float ' + f + ' does not contain a profile '
                         'dir inside its directory. No data will be '
                         'downloaded for this float')
            connection.cwd('..')


        for f in dead_floats:
            to_be_updated = False
            # Update the xml with the current status of the float
            # Check if it must be updated
            f_in_xml = root.findall('wmo_' + str(f))
            if len(f_in_xml) == 0:
                # If this float is new, then add it to the archive
                # and it will be updated
                to_be_updated = True
                f_node = xml_tree.SubElement(root, 'wmo_' + str(f))
                f_node.set('status', 'D')
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
                    connection.cwd(f)
                except:
                    log.info('No directory associated with wmo ' + str(f) +
                             '. This wmo will be skipped!')
                    continue

                _, float_dirs, _ = list_files(connection)
                # Now I look for the profiles dir. This is the folder
                # where all the data are stored
                if 'profiles' in float_dirs:
                    connection.cwd('profiles')
                    float_files, _, perms = list_files(connection)
                    min_len = len(f) + 2
                    to_be_downloaded = [ff for ff in float_files
                                        if len(ff)>min_len
                                        and ff[:min_len] in ['SR'+f,'SD'+f]]
                    if len(to_be_downloaded) > 0:
                        # Copy all file in a local dir with the same name
                        # skipping the one that we already have
                        float_local_dir = join(path, f)
                        ensure_dir(float_local_dir, log, expected = False)
                        for ff in to_be_downloaded:
                            d = download_file(connection, ff, float_local_dir,
                                              log, perms, True)
                            # If the file was downloaded without any problem,
                            # add it to the list of downloaded files
                            if d:
                                downloaded.append(ff)
                    connection.cwd('..')
                else:
                    log.info('Float ' + f + ' does not contain a profile '
                             'dir inside its directory. No data will be '
                             'downloaded for this float')
                connection.cwd('..')
        connection.quit()

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
        that float that starts with 'SR' or 'SD'. Then create a xml file with the
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
        downloaded = []
        A=self.wmo_file_reader()
        wmo_list = A['wmo']
        nFloats=A.size
        wmo_status = dict()
        for l in range(nFloats):
            wmo_status[A[l]['wmo']] = A[l]['status']

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
        connection = FTP(ftp_url)
        connection.login()

        # Enter in the directory tree
        connection.cwd('ifremer/argo/dac/coriolis')

        # Download data for every active float
        for f in wmo_list:

            # Update the xml with the current status of the float
            f_in_xml = root.findall('wmo_' + str(f))
            if len(f_in_xml) == 0:
                f_node = xml_tree.SubElement(root, 'wmo_' + str(f))
                f_node.set('status', wmo_status[f])
            else:
                f_node = [fn for fn in f_in_xml if fn.tag=='wmo_'+str(f)][0]
                f_node.set('status', wmo_status[f])

            try:
                connection.cwd(f)
            except:
                log.info('No directory associated with wmo ' + str(f) +
                         '. This wmo will be skipped!')
                continue

            _, float_dirs, _ = list_files(connection)
            # Now I look for the profiles dir. This is the folder
            # where all the data are stored
            if 'profiles' in float_dirs:
                connection.cwd('profiles')
                float_files, _, perms = list_files(connection)
                min_len = len(f) + 2
                to_be_downloaded = [ff for ff in float_files
                                    if len(ff)>min_len
                                    and ff[:min_len] in ['SR'+f, 'SD'+f]]
                if len(to_be_downloaded) > 0:
                    download_for_f = []
                    # Copy all file in a local dir with the same name
                    float_local_dir = join(path, f)
                    ensure_dir(float_local_dir, log, expected = False)
                    for ff in to_be_downloaded:
                        d = download_file(connection, ff, float_local_dir,
                                          log, perms, skip_if_present)
                        # If the file was downloaded without any problem,
                        # add it to the list of downloaded files
                        if d:
                            downloaded.append(ff)
                            download_for_f.append(ff)
                    if len(download_for_f) == 0:
                        log.info('No updates found for float ' + str(f))                    
                else:
                    log.info('No updates found for float ' + str(f))
                connection.cwd('..')
            else:
                log.info('Float ' + f + ' does not contain a profile '
                         'dir inside its directory. No data will be '
                         'downloaded for this float')
            connection.cwd('..')

        connection.quit()
        
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

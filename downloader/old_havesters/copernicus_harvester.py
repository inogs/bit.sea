# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
from __future__ import print_function

from ftplib import FTP
from os.path import join
from os import listdir

from harvester_interface import HarvesterInterface
from utilities.ftp_utilities import list_files, download_file
from utilities.files_and_dirs import ensure_dir

ftp_url = 'nrt.cmems-du.eu'

user = 'MED_OGS_TRIESTE_IT'
password = 'NEdifupa'

relative_path = "COPERNICUS"


class CopernicusHarvester(HarvesterInterface):
    """
    This is the harvester in charge of download all the files from the
    ftp server medinsitu.hcmr.gr.
    """

    def harvest(self, db_path, log):
        """
        Download all the files inside the remote directories "vessel" and
        "mooring" of the remote ftp server whose modification date is after
        the modification date of the last file in the local dir. Please do
        not put any file in the local directory because this may change the
        date of the last edited file

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

        # Check if the directory for this harvester is present
        # in the database
        path = join(db_path,relative_path)
        ensure_dir(path, log, expected=True)
        # Check if exists the folder "vessel"
        path_vessel = join(path, "vessel")
        ensure_dir(path_vessel, log, expected=True)
        # Check if exists the folder "mooring"
        path_mooring = join(path, "mooring")
        ensure_dir(path_mooring, log, expected=True)

        # Open the connection with the remote archive
        connection = FTP(ftp_url)
        connection.login(user=user, passwd=password)

        # Enter in the folders
        connection.cwd('Core')
        connection.cwd('INSITU_MED_NRT_OBSERVATIONS_013_035')
        connection.cwd('monthly')
        
        # Now I will download everything from the vessel dir
        connection.cwd('vessel')
        log.debug("Entering in dir vessel")

        # Check the last file we have already downloaded
        already_downloaded = listdir(path_vessel)
        file_dates = [int(l.split('_')[1]) for l in already_downloaded]
        if len(file_dates) == 0:
            last_downloaded = 0
        else:
            last_downloaded = max(file_dates)
        log.debug("Last downloaded file on ??/{0:0>2}/{1:0>4}".format(
                   last_downloaded%100, last_downloaded//100))
                       
        # List all the dirs and take only the one that are generated
        # after the last file downloaded
        _, subdirs, _ = list_files(connection)
        subdirs_to_check = [d for d in subdirs if int(d) >= last_downloaded]

        # Download all the file in that dirs
        for d in sorted(subdirs_to_check):
            log.debug("Entering in dir vessel/" + d)
            connection.cwd(d)
            files, _, perms = list_files(connection)
            for f in files:
                if f[:2] == "MO" and f[-3:]==".nc":
                    d = download_file(connection, f, path_vessel,
                                      log, perms, True, False)
                    if d:
                        downloaded.append(f)
            connection.cwd('..')
        connection.cwd('..')

        # Now the same for the mooring dir
        connection.cwd('mooring')
        log.debug("Entering in dir mooring")

        already_downloaded = listdir(path_mooring)
        file_dates = [int(l.split('_')[1]) for l in already_downloaded if l!='incomplete_download.tmp']
        if len(file_dates) == 0:
            last_downloaded = 0
        else:
            last_downloaded = max(file_dates)
        log.debug("Last downloaded file on ??/{0:0>2}/{1:0>4}".format(
                   last_downloaded%100, last_downloaded//100))
        
        _, subdirs, _ = list_files(connection)
        subdirs_to_check = [d for d in subdirs if int(d) >= last_downloaded]

        for d in sorted(subdirs_to_check):
            log.debug("Entering in dir mooring/" + d)
            connection.cwd(d)
            files, _, perms = list_files(connection)
            for f in files:
                if f[:2] == "MO" and f[-3:]==".nc":
                    d = download_file(connection, f, path_mooring,
                                      log, perms, True, False)
                    if d:
                        downloaded.append(f)
            connection.cwd('..')
        connection.cwd('..')
        
        # At the end, download the index
        connection.cwd('..')
        _, _, perms = list_files(connection)
        download_file(connection, 'index_monthly.txt', path,
                      log, perms, False)

        connection.quit()
        return downloaded


    def rebuild(self, db_path, log):
        """
        Download all the files inside the remote directories "vessel" and
        "mooring" of the remote ftp server. If a file already exists, it
        will be rewritten.

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

        # Check if the directory for this harvester is present
        # in the database
        path = join(db_path,relative_path)
        ensure_dir(path, log, expected=False)
        # Check if exists the folder "vessel"
        path_vessel = join(path, "vessel")
        ensure_dir(path_vessel, log, expected=False)
        # Check if exists the folder "mooring"
        path_mooring = join(path, "mooring")
        ensure_dir(path_mooring, log, expected=False)


        # Open the connection with the remote archive
        connection = FTP(ftp_url)
        connection.login(user=user, passwd=password)

        connection.cwd('Core')
        connection.cwd('INSITU_MED_NRT_OBSERVATIONS_013_035')
        connection.cwd('monthly')
        
        # Enter in the folder "vessel"
        connection.cwd('vessel')
        log.debug("Entering in dir vessel")

        # For every subdir, download every netcdf file whose
        # name starts with "MO" and put it in the vessel
        _, subdirs, _ = list_files(connection)
        for d in sorted(subdirs):
            log.debug("Entering in dir vessel/" + d)
            connection.cwd(d)
            files, _, perms = list_files(connection)
            for f in files:
                if f[:2] == "MO" and f[-3:]==".nc":
                    d = download_file(connection, f, path_vessel,
                                      log, perms, False)
                    if d:
                        downloaded.append(f)
            connection.cwd('..')
        connection.cwd('..')
        
        # The same for the other dir
        connection.cwd('mooring')
        log.debug("Entering in dir mooring")

        _, subdirs, _ = list_files(connection)
        for d in sorted(subdirs):
            log.debug("Entering in dir mooring/" + d)
            connection.cwd(d)
            files, _, perms = list_files(connection)
            for f in files:
                if f[:2] == "MO" and f[-3:]==".nc":
                    d = download_file(connection, f, path_mooring,
                                      log, perms, False)
                    if d:
                        downloaded.append(f)
            connection.cwd('..')
        connection.cwd('..')

        # At the end, download the index
        connection.cwd('..')
        _, _, perms = list_files(connection)
        download_file(connection, 'index_monthly.txt', path,
                      log, perms, False)

        connection.quit()
        return downloaded


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

relative_path = "SAT/MULTISENSOR/1Km/NRT/DAILY/ORIG/"


class Sat_ms1kmNRT_Harvester(HarvesterInterface):
    """
    This is the harvester in charge of download all the files from the
    ftp server myocean.artov.isac.cnr.it.
    """
    def harvest(self, db_path, log):
        """
        Download all the files inside a remote directory of the ftp server
        whose modification date is after the modification date of the last
        file in the local dir. Download the files if the contain
        '-NRT-' in their name.

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

        # Open the connection with the remote archive
        connection = FTP(ftp_url)
        connection.login(user=user, passwd=password)


        connection.cwd('Core')
        connection.cwd('OCEANCOLOUR_MED_CHL_L3_NRT_OBSERVATIONS_009_040')
        connection.cwd('dataset-oc-med-chl-multi-l3-chl_1km_daily-rt-v02')
        
        # List all the local files
        loc_files = [f for f in listdir(path) if f !='incomplete_download.tmp']
        
        # If there are no files, download everything
        if len(loc_files)==0:
            log.info('No local files found! Everything will be '
                     'downloaded from the remote repository!')
            _, years, _ = list_files(connection)
            for year in years:
                connection.cwd(year)
                _, months, _ = list_files(connection)
                for month in months:
                    connection.cwd(month)
                    files, _, perms = list_files(connection)
                    files_to_be_downloaded = [f for f in files if '-NRT-' in f]
                    for f in files_to_be_downloaded:
                        d = download_file(connection, f, path,
                                               log, perms, False)
                        if d:
                            downloaded.append(f)
                    connection.cwd('..')
                connection.cwd('..')
        else:
            loc_files.sort()
            last_file = loc_files[-1]
            last_year = int(last_file[0:4])
            _, years, _ = list_files(connection)
            new_years = [y for y in years if int(y)>last_year]
            # Enter in the folder with the year of the last downloaded
            # file and download every file which is newer than that
            connection.cwd(str(last_year))
            _, months, _ = list_files(connection)
            for month in months:
                connection.cwd(month)
                files, _, perms = list_files(connection)
                files_to_be_downloaded = [f for f in files if '-NRT-' in f]
                for f in files_to_be_downloaded:
                    if f > last_file:
                        d = download_file(connection, f, path,
                                          log, perms, True, True)
                        if d:
                            downloaded.append(f)
                connection.cwd('..')
            connection.cwd('..')    
            # Now we will download what is in the folders of the years
            # after the last file
            for year in new_years:
                connection.cwd(year)
                _, months, _ = list_files(connection)
                for month in months:
                    connection.cwd(month)
                    files, _, perms = list_files(connection)
                    files_to_be_downloaded = [f for f in files if '-NRT-' in f]
                    for f in files_to_be_downloaded:
                        d = download_file(connection, f, path,
                                          log, perms, True, True)
                        if d:
                            downloaded.append(f)
                    connection.cwd('..')
            connection.cwd('..')
            # Warning if we found a lot of updates or no updates at all
            if len(downloaded) == 0:
                log.info('No updates found!')
            if len(downloaded) >1 : 
                warn_message = 'Downloaded more than one file:'
                for f in downloaded:
                    warn_message += '\n   - ' + str(f)
                log.info(warn_message, split_lines=False)

        connection.quit()
        return downloaded


    def rebuild(self, db_path, log):
        """
        Download all the files inside a remote directory of the ftp server. If the
        file is already present on the local directory, rewrite it. Do not download
        the file that contain the string '-NRT-' in their filename.

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

        # Open the connection with the remote archive
        connection = FTP(ftp_url)
        connection.login(user=user, passwd=password)

        connection.cwd('Core')
        connection.cwd('OCEANCOLOUR_MED_CHL_L3_NRT_OBSERVATIONS_009_040')
        connection.cwd('dataset-oc-med-chl-multi-l3-chl_1km_daily-rt-v02')
        
        _, years, _ = list_files(connection)
        for year in years:
            connection.cwd(year)
            _, months, _ = list_files(connection)
            for month in months:
                connection.cwd(month)
                files, _, perms = list_files(connection)
                files_to_be_downloaded = [f for f in files if '-NRT-' in f]
                for f in files_to_be_downloaded:
                    d = download_file(connection, f, path,
                                           log, perms, False)
                    if d:
                        downloaded.append(f)
                connection.cwd('..')
            connection.cwd('..')
        connection.quit()
        return downloaded

from __future__ import print_function

from ftplib import FTP
from os.path import join, exists
from os import listdir

from utilities.ftp_utilities import list_files
from utilities.files_and_dirs import ensure_dir


ftp_url = 'myocean.artov.isac.cnr.it'

user = 'MED_OGS_TRIESTE_IT'
password = 'NEdifupa'
relative_path = "myocean"

class UnexpectedDir(Exception):
    pass

class MyOceanHarvester(object):

    def harvest(self, db_path, log):
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

        # Enter in the folder "Intermediate"
        connection.cwd('Intermediate')

        # Enter in "OCEANCOLOUR_MED_CHL_L4_NRT_OBSERVATIONS_009_060"
        connection.cwd('OCEANCOLOUR_MED_CHL_L4_NRT_OBSERVATIONS_009_060')

        # Enter in "dataset-oc-med-chl-modis_a-l4-chl_7km_daily-rt-v02"
        connection.cwd('dataset-oc-med-chl-modis_a-l4-chl_7km_daily-rt-v02')
        
        # List all the local files
        loc_files = listdir(path)
        
        # If there are no files, download everything
        if len(loc_files)==0:
            log.queerness('No local files found! Everything will be '
                          ' downloaded from the remote repository!')
            _, years, _ = list_files(connection)
            for year in years:
                connection.cwd(year)
                files, _, perms = list_files(connection)
                for f in files:
                    d = self.download_file(connection, f, path,
                                           perms, log, False)
                    if d:
                        downloaded.append(f)
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
            files, _, perms = list_files(connection)
            for f in files:
                if f > last_file:
                    d = self.download_file(connection, f, path,
                                           perms, log, True, True)
                    if d:
                        downloaded.append(f)
            connection.cwd('..')
            # Now we will download what is in the folders of the years
            # after the last file
            for year in new_years:
                connection.cwd(year)
                files, _, perms = list_files(connection)
                for f in files:
                    d = self.download_file(connection, f, path,
                                           perms, log, True, True)
                    if d:
                        downloaded.append(f)
                connection.cwd('..')

            # Warning if we found a lot of updates or no updates at all
            if len(downloaded) == 0:
                log.queerness('No updates found!')
            if len(downloaded) >1 : 
                warn_message = 'Downloaded more than one file:'
                for f in downloaded:
                    warn_message += '\n   - ' + str(f)
                log.queerness(warn_message, split_lines=False)

        return downloaded


    def rebuild(self, db_path, log):
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

        # Enter in the folder "Intermediate"
        connection.cwd('Intermediate')

        # Enter in "OCEANCOLOUR_MED_CHL_L4_NRT_OBSERVATIONS_009_060"
        connection.cwd('OCEANCOLOUR_MED_CHL_L4_NRT_OBSERVATIONS_009_060')

        # Enter in "dataset-oc-med-chl-modis_a-l4-chl_7km_daily-rt-v02"
        connection.cwd('dataset-oc-med-chl-modis_a-l4-chl_7km_daily-rt-v02')
        
        _, years, _ = list_files(connection)
        for year in years:
            connection.cwd(year)
            files, _, perms = list_files(connection)
            for f in files:
                d = self.download_file(connection, f, path,
                                       perms, log, False)
                if d:
                    downloaded.append(f)
            connection.cwd('..')
        return downloaded




    def download_file(self, connection, f, path, perms, log,
                      skip_if_exists=True, skip_is_strange = False):
        if skip_if_exists:
            # Check if the file already exists
            if exists(join(path, f)):
                if skip_is_strange:
                    log.queerness('Skipping file ' + f + ' because'
                                 ' it was already downloaded!')
                else:
                    log.info('Skipping file ' + f + ' because '
                             'it was already downloaded')
                return False

        # Check if the file is readable
        if 'read' not in perms[f]['others']:
            log.queerness('The file ' + f + ' is not '
                          'readable. Check the '
                          'permissions')
            return False

        # Download the file
        log.info('Downloading ' + f + '...')
        with open(join(path, f), 'w') as loc_file:
            connection.retrbinary('RETR ' + str(f), loc_file.write)
        log.info('File ' + f + ' downloaded!')
        return True


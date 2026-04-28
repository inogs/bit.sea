import argparse
import os
from ftplib import FTP
from pathlib import Path, PurePosixPath
def argument():
    parser = argparse.ArgumentParser(description = '''
    Downloads:
    - netcdf files of argo floats in Med_floats.txt if not already present
    - netcdf files of update file, in force mode
                                 
    The script also removes files that begin with "SR" when "SD" files with the same WMO number and rise up number exist.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--remotedir',"-r",
                                type = str,
                                required = True,
                                help = 'remote directory of the argo indexed files',
                                default = '/ifremer/argo/dac/')
    parser.add_argument(   '--coriolisdir',"-c",
                                type = str,
                                required = True,
                                help = 'local directory of the argo indexed files',
                                default = '/leonardo_work/OGS_prod2528_0/OPA/V12C/ONLINE/CORIOLIS/')
    parser.add_argument( '--update_file',"-u",
                                type = str,
                                required = True,
                                help = 'path to the update file containing the list of floats to download',
                                default = '/leonardo_work/OGS_prod2528_0/OPA/V12C/ONLINE/CORIOLIS/download/DIFF_floats.txt')
    parser.add_argument(   '--indexer_file',"-i",
                                type = str,
                                required = False,
                                help = 'path to the indexer file containing the list of floats that have already been downloaded',
                                default = '/leonardo_work/OGS_prod2528_0/OPA/V12C/ONLINE/CORIOLIS/download/Med_floats.txt')

    return parser.parse_args()

args = argument()

remote_dir = PurePosixPath(args.remotedir)
local_coriolis_dir = Path(args.coriolisdir)
update_file = args.update_file
if args.indexer_file:
    indexer_file = args.indexer_file
else:
    indexer_file = None

local_tmp_download_path = local_coriolis_dir / "download" / "tmp"

def parse_float_line(line):
    split_line = line.split(",")

    if len(split_line) > 0:
        float_path_remote = PurePosixPath(split_line[0])
        float_path_remote_split = float_path_remote.parts
        wmo = float_path_remote_split[1]
        float_file_name = float_path_remote_split[-1]
        float_path_local = local_coriolis_dir / wmo / float_file_name
        return float_path_remote, float_path_local, float_file_name
    else:
        return None, None, None

def download_move_float(float_path_remote, float_path_local, float_file_name):
    local_tmp_download_path.mkdir(parents=True, exist_ok=True)
    float_path_local.parent.mkdir(parents=True, exist_ok=True)
    tmp_file_path = local_tmp_download_path / float_file_name
    remote_file_path = str(remote_dir / float_path_remote)

    print(f"Downloading {remote_file_path} from ftp.ifremer.fr")

    try:
        with FTP("ftp.ifremer.fr") as connection:
            connection.login()
            with open(tmp_file_path, "wb") as tmp_file:
                connection.retrbinary(f"RETR {remote_file_path}", tmp_file.write)
    except Exception as exc:
        if tmp_file_path.exists():
            os.remove(tmp_file_path)
        print(f"Failed to download {float_file_name}: {exc}")
        return False

    os.rename(tmp_file_path, float_path_local)
    print(f"Successfully downloaded and moved {float_file_name} to {float_path_local}")
    return True

# download from indexer file list (if provided)
if indexer_file is not None:
    with open(indexer_file, "r") as f:
        already_downloaded_floats = set(line.strip() for line in f)
    for line in already_downloaded_floats:
        float_path_remote, float_path_local, float_file_name = parse_float_line(line)
        if float_path_local is not None and \
            float_path_remote is not None and \
                float_file_name is not None and \
                    not os.path.exists(float_path_local):
            download_move_float(float_path_remote, float_path_local, float_file_name)

# download from DIFF_floats.txt list
with open(update_file, "r") as f:
    floats_to_download = set(line.strip() for line in f)
for line in floats_to_download:
    float_path_remote, float_path_local, float_file_name = parse_float_line(line)
    if float_path_local is not None and \
        float_path_remote is not None and \
            float_file_name is not None:
        download_move_float(float_path_remote, float_path_local, float_file_name)

# clean up local Coriolis directory by removing files 
# that begin with "SR" when "SD" files with the same 
# WMO number and rise up number exist.
if os.path.isdir(local_coriolis_dir):
    for wmo_dir in os.listdir(local_coriolis_dir):
        local_float_dir = local_coriolis_dir / wmo_dir
        if (not local_float_dir.is_dir()) or (len(wmo_dir) != 7) or (not wmo_dir.isdigit()):
            continue
        for filename in os.listdir(local_float_dir):
            if filename.startswith('SD') and filename.endswith('.nc'):
                sr_path = local_float_dir / ('SR' + filename[2:])
                if sr_path.exists():
                    os.remove(sr_path)
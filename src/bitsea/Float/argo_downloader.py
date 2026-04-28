import argparse
import os
import subprocess as sp
from pathlib import Path
def argument():
    parser = argparse.ArgumentParser(description = '''
    Downloads:
    - netcdf files of argo floats in Med_floats.txt if not already present
    - netcdf files of update file, in force mode
                                 
    The script also removes files that begin with "SR" when "SD" files with the same WMO number and rise up number exist.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--getcommand',"-g",
                                type = str,
                                required = False,
                                help = 'path of the ncftpget command',
                                default = '/leonardo/home/usera07ogs/a07ogs00/OPA/V12C-dev/HOST/leonardo/bin/ncftpget')
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

ncftpget_command = args.getcommand
remote_dir = Path("/ifremer/argo/dac/")
local_coriolis_dir = Path(args.coriolisdir)
update_file = args.update_file
if args.indexer_file:
    indexer_file = args.indexer_file
else:
    indexer_file = None

local_tmp_download_path = local_coriolis_dir / "download" / "tmp"

def parse_float_line(line):
    """
    Parses one line of the argo floats list.
    Input: 
    - line : str
    Output:
    - float_path_remote : str or None  remote path (source file to be downloaded)
    - float_path_local  : str or None  local path  (destination)
    - float_file_name   : str or None  file name
    """
    split_line = line.split(",")

    if len(split_line) > 0:
        float_path_remote = Path(split_line[0])
        float_path_remote_split = float_path_remote.parts
        float_dir = float_path_remote_split[1]   #wmo number
        float_file_name = float_path_remote_split[-1]
        float_path_local = local_coriolis_dir / float_dir / float_file_name
        return float_path_remote, float_path_local, float_file_name
    else:
        return None, None, None
    
def download_move_float(float_path_remote, float_path_local, float_file_name):
    """
    Downloads the netcdf file containing the profiles associated to a specified argo float and rise up number.
    Moves the downloaded file to the specified local directory.
    Input:
    - float_path_remote : str or None  remote path (source file to be downloaded)
    - float_path_local  : str or None  local path  (destination)
    - float_file_name   : str or None  file name
    Output:
    - bool              : False/True   failure/success output states
    """
    remote_file_path = str(remote_dir / float_path_remote)
    command = f"{ncftpget_command} -u anonymous ftp.ifremer.fr {local_tmp_download_path} {remote_file_path)}"
    
    print(f"Downloading {remote_file_path} from ftp.ifremer.fr")

    res=sp.run(command, shell=True, stdout=sp.PIPE, stderr=sp.STDOUT)
    res_str = res.stdout.decode("utf-8")
    if res.returncode == 0:
        os.rename(local_tmp_download_path + "/" + float_file_name, float_path_local)
        print(f"Download OK: {float_path_local}")
    else:
        print(f"""Download FAILED: {float_file_name}.
              Command output: {res_str}""")
    return res.returncode == 0


# download from indexer file list (if provided)
print("Downloading missing argo float profiles (sync list with existing files)...")
if indexer_file is not None:
    with open(indexer_file, "r") as f:
        already_downloaded_floats = set(line.strip() for line in f)
    for line in already_downloaded_floats:
        float_path_remote, float_path_local, float_file_name = parse_float_line(line)
        parser_success = float_path_local is not None and float_path_remote is not None and float_file_name is not None
        if float_path_local is not None:
            profile_is_missing = not os.path.exists(float_path_local)
        else:
            profile_is_missing = False
        if parser_success and profile_is_missing:
            download_move_float(float_path_remote, float_path_local, float_file_name)
print("...Done.")

# download from DIFF_floats.txt list
print("Downloading new argo float profiles...")
with open(update_file, "r") as f:
    floats_to_download = set(line.strip() for line in f)
for line in floats_to_download:
    float_path_remote, float_path_local, float_file_name = parse_float_line(line)
    parser_success = float_path_local is not None and float_path_remote is not None and float_file_name is not None
    if parser_success:
        download_move_float(float_path_remote, float_path_local, float_file_name)
print("...Done.")

# clean up local Coriolis directory by removing files 
# that begin with "SR" when "SD" files with the same 
# WMO number and rise up number exist.
print("Looking for and removing useless SR argo profiles...")
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
                    print(f"REMOVED: {str(sr_path)}")
print("...Done.")
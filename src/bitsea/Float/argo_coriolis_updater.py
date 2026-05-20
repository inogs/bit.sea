import argparse
import gzip
import os
import shutil
from datetime import datetime
from ftplib import FTP
from pathlib import Path

import numpy as np

from bitsea.basins import V2
from bitsea.commons.utils import file2stringlist

def argument():
    parser = argparse.ArgumentParser(description = '''
    UPDATE SCRIPT FOR ARGO FLOATS IN CORIOLIS

    - removes the old argo_synthetic-profile_index.txt.gz if exists
    - downloads argo_synthetic-profile_index.txt.gz from ifremer ftp server
    - extracts the txt file from the gz file
                                     
    - netcdf files of argo floats in Med_floats.txt if not already present
    - netcdf files of update file, in force mode
                                 
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--coriolis_dl_dir',"-c",
                                type = str,
                                required = True,
                                help = 'local directory of the argo indexed files',
                                default = '/leonardo_work/OGS_prod2528_0/OPA/V12C/ONLINE/CORIOLIS/download/')
    parser.add_argument( '--update_file',"-u",
                                type = str,
                                required = True,
                                help = 'path to the update file containing the list of floats to download',
                                default = '/leonardo_work/OGS_prod2528_0/OPA/V12C/ONLINE/CORIOLIS/download/DIFF_floats.txt')
    parser.add_argument(   '--indexer_file',"-i",
                                type = str,
                                required = True,
                                help = 'path to the indexer file containing the list of floats that have already been downloaded',
                                default = '/leonardo_work/OGS_prod2528_0/OPA/V12C/ONLINE/CORIOLIS/download/Med_floats.txt')

    return parser.parse_args()

args = argument()


def build_mediterranean_index(input_file: Path | str, output_file: Path | str) -> None:
    input_file = str(input_file)
    output_file = str(output_file)

    input_dtype = np.dtype([
        ('file_name', 'S200'),
        ('date', 'S200'),
        ('latitude', np.float32),
        ('longitude', np.float32),
        ('ocean', 'S10'),
        ('profiler_type', int),
        ('institution', 'S10'),
        ('parameters', 'S200'),
        ('parameter_data_mode', 'S100'),
        ('date_update', 'S200')
    ])

    index_file = np.loadtxt(input_file, dtype=input_dtype, delimiter=",", ndmin=1, skiprows=9)

    # file,date,latitude,longitude,ocean,profiler_type,institution,parameters,parameter_data_mode,date_update
    name_file = index_file['file_name']
    list_name = []
    for name in name_file:
        list_name.append(name.decode().startswith('coriolis'))

    coriolis = index_file[list_name]

    mediterr = V2.med

    lat_lon_list = []
    for ele in coriolis:
        lat_lon_list.append(mediterr.is_inside(lat=ele['latitude'], lon=ele['longitude']))

    output_dtype = np.dtype([
        ('file_name', 'U200'),
        ('date', 'U200'),
        ('latitude', np.float32),
        ('longitude', np.float32),
        ('ocean', 'U10'),
        ('profiler_type', int),
        ('institution', 'U10'),
        ('parameters', 'U200'),
        ('parameter_data_mode', 'U100'),
        ('date_update', 'U200')
    ])

    n_float_med = sum(lat_lon_list)
    med_float = np.zeros((n_float_med,), dtype=output_dtype)
    med_float[:] = coriolis[lat_lon_list]

    np.savetxt(output_file, med_float, fmt="%s,%s,%f,%f,%s,%d,%s,%s,%s,%s")

def build_difference_file(input_new: Path | str, input_old: Path | str, output_file: Path | str) -> None:
    old_list = file2stringlist(str(input_old))
    new_list = file2stringlist(str(input_new))

    diff = []
    for line in new_list:
        if line not in old_list:
            diff.append(line)

    lines = []
    for line in diff:
        lines.append(line + "\n")

    with open(output_file, "wt") as file_handle:
        file_handle.writelines(lines)


def download_argo_index_file(c_dl_dir:Path):
    """
    Downloads the argo synthetic profile index file.
    Input:
    - c_dl_dir : Path  local directory where the file will be downloaded
    Output:
    - bool              : False/True   failure/success output states
    """

    file_name = "argo_synthetic-profile_index.txt.gz"
    local_file_path = c_dl_dir / file_name
    remote_file_path = "/ifremer/argo/" + file_name

    print(f"Downloading {remote_file_path} from ftp.ifremer.fr")

    try:
        with FTP("ftp.ifremer.fr") as connection:
            connection.login()
            with open(local_file_path, "wb") as file:
                connection.retrbinary(f"RETR {remote_file_path}", file.write)
    except Exception as exc:
        if local_file_path.exists():
            os.remove(local_file_path)
        print(f"""Download FAILED: {file_name}.
              Command output: {exc}""")
        return False
    return True
    
c_dl_dir = Path(args.coriolis_dl_dir)
update_file = Path(args.update_file) if args.update_file else None
indexer_file = Path(args.indexer_file) if args.indexer_file else None

if not c_dl_dir.exists():
    print(f"Creating directory {c_dl_dir}")
    c_dl_dir.mkdir(parents=True, exist_ok=True)
if update_file is None or indexer_file is None:
    print("Update file or indexer file not provided. Exiting.")
    exit(1)

gz_file = c_dl_dir / "argo_synthetic-profile_index.txt.gz"
txt_file = c_dl_dir / "argo_synthetic-profile_index.txt"

if gz_file.exists():
    os.remove(gz_file)

if not download_argo_index_file(c_dl_dir):
    print("Cannot continue: failed to download argo index file.")
    raise SystemExit(1)

with gzip.open(gz_file, "rb") as f_in, open(txt_file, "wb") as f_out:
    shutil.copyfileobj(f_in, f_out)


# correct file from blank parts
corr_txt_file = c_dl_dir / "argo_synthetic-profile_index_CORR.txt"
with open(txt_file, "r", encoding="utf-8", errors="ignore") as src, open(corr_txt_file, "w", encoding="utf-8") as dst:
    for line in src:
        if ",," not in line and "D.nc" not in line:
            dst.write(line)

#launch the reader -> obtain the mediterranean floater file
build_mediterranean_index(corr_txt_file, indexer_file)
#matching old vs new -> return a file called DIFF_floats.txt in which there is the list of floats that are different between the two files
old_indexer_file = c_dl_dir / "Med_floats_OLD.txt"
build_difference_file(indexer_file, old_indexer_file, update_file)

daily_log_dir = Path(os.environ["OPA_LOGDIR"]) / "daily"
if not daily_log_dir.exists():
    print(f"Creating directory {daily_log_dir}")
    daily_log_dir.mkdir(parents=True, exist_ok=True)

diff_log_file = daily_log_dir / f"DIFF_floats.{datetime.now().strftime('%Y%m%d-%H:%M:%S')}.txt"
shutil.copyfile(update_file, diff_log_file)



import argparse
from pathlib import Path

import numpy as np

from bitsea.basins import V2


def argument():
    parser = argparse.ArgumentParser(description = '''
    Read file argo, gives back file with mediterannean argo floats
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--inputfile',"-i",
                                type = str,
                                required = True,
                                help = 'input file argo txt, eg=')
    parser.add_argument(   '--outfile',"-o",
                                type = str,
                                required = True,
                                help = 'output file argo')
    return parser.parse_args()


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


def build_mediterranean_index_from_dir(input_corr_index_file: Path | str, indexer_file: Path | str) -> None:

    build_mediterranean_index(input_file=input_corr_index_file, output_file=indexer_file)


def main() -> None:
    args = argument()
    build_mediterranean_index(input_file=args.inputfile, output_file=args.outfile)


if __name__ == "__main__":
    main()

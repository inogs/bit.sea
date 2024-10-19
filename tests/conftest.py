from os import getenv
from pathlib import Path
from tempfile import TemporaryDirectory
from tempfile import NamedTemporaryFile
from typing import Optional
from warnings import warn
import pytest
import requests
import shutil
import tarfile

from bitsea.utilities.argparse_types import existing_dir_path
from bitsea.utilities.argparse_types import existing_file_path

TEST_DATA_DIR_VAR = 'BITSEA_TEST_DATA_DIR'
TEST_DATA_FILE_VAR = 'BITSEA_TEST_DATA_FILE'
TEST_DATA_URL = "https://medeaf.ogs.it/internal-validation/spiani00/bitsea_test_data.tar.xz"


def pytest_addoption(parser):
    parser.addoption(
        "--test-data-dir",
        required=False,
        type=existing_dir_path,
        default=None
    )
    parser.addoption(
        "--test-data-file",
        required=False,
        type=existing_file_path,
        default=None
    )


def decompress_data_file_internal(data_file: Path):
    output_dir = TemporaryDirectory()

    output_dir_path = Path(output_dir.name)
    print(f'Decompressing file {data_file} inside dir {output_dir_path}...')
    with tarfile.open(data_file) as f:
        f.extractall(output_dir_path)
    print('Done!')

    return output_dir


def decompress_data_file(data_file: Optional[Path]):
    if data_file is None:
        with NamedTemporaryFile(mode="w+b", suffix='.tar.xz') as f:
            f_path = Path(f.name)
            print(f'Downloading from {TEST_DATA_URL} into file {f_path}')
            if data_file is None:
                r = requests.get(TEST_DATA_URL, stream=True)
                for data in r.iter_content(chunk_size=1024 * 512):
                    f.write(data)
            print('Done!')
            return decompress_data_file_internal(f_path)
    return decompress_data_file_internal(data_file)


def pytest_sessionstart(session):
    delete_test_data_dir = False

    # We try to find the test data dir from the command line
    test_data_dir = session.config.option.test_data_dir
    test_data_file = session.config.option.test_data_file

    if test_data_dir is not None and test_data_file is not None:
        warn(
            '--test-data-dir and --test-data-file are mutually exclusive. The '
            'second will be ignored!'
        )

    # If it has not been submitted, we try from the env variables
    if test_data_dir is None:
        if env_var_content := getenv(TEST_DATA_DIR_VAR) is not None:
            test_data_dir = Path(env_var_content)

    if test_data_dir is not None and test_data_file is not None:
        warn(
            f'The submitted value for --test-data-file will be ignored because '
            f'the directory with the test data is being obtained from the ENV '
            f'variable {TEST_DATA_DIR_VAR}'
        )

    # If it is still None, we build it, but then we must remember to delete it
    # after the execution of the tests
    if test_data_dir is None:
        delete_test_data_dir = False

        if test_data_file is None:
            if env_var_content := getenv(TEST_DATA_FILE_VAR) is not None:
                test_data_file = Path(env_var_content)

        test_data_dir = decompress_data_file(test_data_file)

    session.config.option.test_data_dir = test_data_dir
    session.config.option.delete_test_data_dir = delete_test_data_dir


def pytest_sessionfinish(session):
    if session.config.option.delete_test_data_dir:
        try:
            session.config.option.test_data_dir.cleanup()
        except AttributeError:
            shutil.rmtree(session.config.option.test_data_dir)


@pytest.fixture
def test_data_dir(request) -> Path:
    try:
        return Path(request.session.config.option.test_data_dir.name)
    except AttributeError:
        return request.session.config.option.test_data_dir

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

# Configuration for finding the data for the tests;
# Specifies the name of the environment variable that describes what
# if the path of the directory with the data for the tests
TEST_DATA_DIR_VAR = 'BITSEA_TEST_DATA_DIR'

# If the directory does not exist, there is another environment
# variable that specifies where to find a compress file with the
# data for the tests
TEST_DATA_FILE_VAR = 'BITSEA_TEST_DATA_FILE'

# If also the file does not exist, we can still download it from
# the following URL
TEST_DATA_URL = "https://medeaf.ogs.it/internal-validation/spiani00/bitsea_test_data.tar.xz"

# If this directory exists, then we do not download the file from the
# previous URL, but we use this directory instead
DEFAULT_TEST_DATA_DIR = Path(__file__).parent / 'data'


def pytest_addoption(parser):
    # Add two arguments to the cli of pytest to specify
    parser.addoption(
        "--test-data-dir",
        required=False,
        type=existing_dir_path,
        default=None,
        help=f"The path of the directory with the files for the tests. "
             f"If this argument is not submitted, the code will use a "
             f"temporary directory and will generate the files by "
             f"extracting a tar.xz file inside this directory before "
             f"running the tests. This is equivalent to set a environment "
             f"variable named {TEST_DATA_DIR_VAR}; if this argument is "
             f"submitted while the env var {TEST_DATA_DIR_VAR} is defined, "
             f"the value of the env variable will be ignored."
    )
    parser.addoption(
        "--test-data-file",
        required=False,
        type=existing_file_path,
        default=None,
        help=f"The path of the compressed file that must be "
             f"decompressed to extract the data for the tests. This "
             f"argument is ignored if the code is aware of a path to a "
             f"directory where the data is already stored (either "
             f"because the --test-data-dir is specified or because the "
             f"env variable {TEST_DATA_DIR_VAR} is defined)."
    )


def decompress_data_file_internal(data_file: Path):
    """
    Decompresses the data_file inside a TemporaryDirectory and returns
    the TemporaryDirectory object as output.
    """
    output_dir = TemporaryDirectory()

    output_dir_path = Path(output_dir.name)
    print(f'Decompressing file {data_file} inside dir {output_dir_path}...')
    with tarfile.open(data_file) as f:
        f.extractall(output_dir_path)
    print('Done!')

    return output_dir


def decompress_data_file(data_file: Optional[Path]):
    """
    Decompresses the data_file inside a TemporaryDirectory and returns
    the TemporaryDirectory object as output. If `data_file` is None, it will
    be downloaded from an external resource, and it will be deleted as soon as
    it has been decompressed.
    """
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


def prepare_test_data_dir(session):
    """
        Prepares the test data directory using a multistep process:

        The function follows this logic:
        1. Check if the test data directory path was provided via the
           command line. If so, use this directory.
        2. If not, check if the directory path was supplied through an
           environment variable. If available, use this directory.
        3. If neither is provided, check if the `DEFAULT_TEST_DATA_DIR`
           exists; if it does, use this directory.
        4. Check if a compressed test data file path was provided via
           the command line. If so, decompress the file to a temporary
           directory and use that directory.
        5. Similarly, check if the compressed file path was provided
           through an environment variable. If found, decompress the
           file to a temporary directory and use it.
        6. If all the above checks fail, download the compressed test
           data file from the external URL specified by
           `TEST_DATA_URL`, decompress it into a temporary directory,
           and use it.

        If a temporary directory is used, it will be deleted after the
        test run.

        This function sets two configuration parameters for the
        session:
        - `test_data_dir`: the path to the test data directory
          (either a standard directory or a `TemporaryDirectory`).
        - `delete_test_data_dir`: a boolean indicating whether
           the test data directory should be deleted after the tests
           are complete.
    """
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
        if (env_var_content := getenv(TEST_DATA_DIR_VAR)) is not None:
            test_data_dir = Path(env_var_content)

    if test_data_dir is not None and test_data_file is not None:
        warn(
            f'The submitted value for --test-data-file will be ignored because '
            f'the directory with the test data is being obtained from the ENV '
            f'variable {TEST_DATA_DIR_VAR}'
        )

    if DEFAULT_TEST_DATA_DIR.is_dir():
        test_data_dir = DEFAULT_TEST_DATA_DIR

    # If it is still None, we build it, but then we must remember to delete it
    # after the execution of the tests
    if test_data_dir is None:
        delete_test_data_dir = False

        if test_data_file is None:
            if (env_var_content := getenv(TEST_DATA_FILE_VAR)) is not None:
                test_data_file = Path(env_var_content)

        test_data_dir = decompress_data_file(test_data_file)

    session.config.option.delete_test_data_dir= delete_test_data_dir
    session.config.option.test_data_dir = test_data_dir


def pytest_sessionstart(session):
    prepare_test_data_dir(session)


def pytest_sessionfinish(session):
    if session.config.option.delete_test_data_dir:
        try:
            session.config.option.test_data_dir.cleanup()
        except AttributeError:
            shutil.rmtree(session.config.option.test_data_dir)


@pytest.fixture
def test_data_dir(request) -> Path:
    if not isinstance(request.session.config.option.test_data_dir, Path):
        return Path(request.session.config.option.test_data_dir.name)
    return request.session.config.option.test_data_dir

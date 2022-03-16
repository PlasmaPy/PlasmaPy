import numpy as np
import os
import pytest

from plasmapy.utils.data import downloader

test_files = [
    # Test downloading a file
    ("NIST_PSTAR_aluminum.txt", None),
    # Test with a different file type
    ("plasmapy_logo.png", None),
    # Test an h5 file
    ("test.h5", None),
    # Test a file without an extension raises an exception
    ("missing_an_extension", ValueError),
    # Test that trying to download a file that doesn't exist raises an
    # exception.
    ("not-a-real-file.txt", OSError),
]


@pytest.mark.parametrize("filename,expected", test_files)
def test_get_file(filename, expected, tmp_path):
    """
    Test the get_file function

    """
    # Delete file if it already exists, so the test always downloads it
    dl_path = os.path.join(tmp_path, filename)
    if os.path.exists(dl_path):
        os.remove(dl_path)

    if expected is not None:
        with pytest.raises(expected):
            path = downloader.get_file(filename, directory=tmp_path)

    else:
        # Download data (or check that it already exists)
        path = downloader.get_file(filename, directory=tmp_path)

        # Get the file again, already existing so it doesn't download it again
        path = downloader.get_file(filename, directory=tmp_path)


def test_get_file_NIST_PSTAR_datafile(tmp_path):
    """
    Test the get_file function on a NIST PSTAR datafile

    """
    filename = "NIST_PSTAR_aluminum.txt"

    # Delete file if it already exists, so the test always downloads it
    dl_path = os.path.join(tmp_path, filename)
    if os.path.exists(dl_path):
        os.remove(dl_path)

    # Download data (or check that it already exists)
    path = downloader.get_file(filename, directory=tmp_path)

    arr = np.loadtxt(path, skiprows=7)
    assert np.allclose(arr[0, :], np.array([1e-3, 1.043e2]))

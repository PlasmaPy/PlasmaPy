from pathlib import Path

import numpy as np
import pytest

from plasmapy.utils.data.downloader import Downloader


@pytest.fixture()
def downloader(tmp_path):
    return Downloader(directory=tmp_path)


test_files = [
    # Test downloading a file
    ("NIST_PSTAR_aluminum.txt", None),
    # Test with a different file type
    ("plasmapy_logo.png", None),
    # Test an h5 file
    ("test.h5", None),
    # Test that trying to download a file that doesn't exist raises an
    # exception.
    ("not-a-real-file.txt", ValueError),
]


@pytest.mark.parametrize(("filename", "expected"), test_files)
def test_get_file(filename, expected, downloader) -> None:
    """Test the get_file function."""

    if expected is not None:
        with pytest.raises(expected):
            downloader.get_file(filename)
    else:
        # Download data (or check that it already exists)
        downloader.get_file(filename)

        # Get the file again, already existing so it doesn't download it again
        downloader.get_file(filename)


def test_get_file_NIST_PSTAR_datafile(downloader) -> None:
    """Test the get_file function on a NIST PSTAR datafile."""
    # Download data (or check that it already exists)
    path = downloader.get_file("NIST_PSTAR_aluminum.txt")

    arr = np.loadtxt(path, skiprows=7)
    assert np.allclose(arr[0, :], np.array([1e-3, 1.043e2]))


test_urls = [
    # Test with a page we know is up if the tests are running
    ("https://github.com/PlasmaPy/PlasmaPy", None),
    # Test with a known 404
    ("https://www.google.com/404", ValueError),
]


@pytest.mark.parametrize(("url", "expected"), test_urls)
def test_http_request(downloader, url, expected):
    """
    Test exceptions from http downloader
    """
    if expected is None:
        downloader._http_request(url)
    else:
        with pytest.raises(expected):
            downloader._http_request(url)


def test_multiple_resource_calls(tmp_path, downloader):
    """
    Test various file retrieval modes
    """
    # Create a dummy file
    filename = "NIST_PSTAR_aluminum.txt"

    # Download from repository
    downloader.get_file(filename)

    # Return local copy, since file now exists
    downloader.get_file(filename)

    # Retrieve a local file that isn't on the remote
    # First create the file
    filename2 = "not_on_the_repo.txt"
    filepath2 = Path(tmp_path, filename2)
    with filepath2.open("w") as f:
        f.write("Not data")
    # Add it to the blob file
    downloader._local_blob_dict[filename2] = "sha"
    downloader._write_blobfile()
    # Now try retrieving it
    downloader.get_file(filename2)

    # Error is raised when a file isn't local or on the remote
    with pytest.raises(ValueError):
        downloader.get_file("not_anywhere.txt")

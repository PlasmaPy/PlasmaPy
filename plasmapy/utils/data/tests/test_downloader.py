from pathlib import Path

import numpy as np
import pytest

from plasmapy.utils.data.downloader import Downloader

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
@pytest.mark.flaky(reruns=5)  # in case of intermittent connection to World Wide Web™
def test_get_file(filename, expected, tmp_path) -> None:
    """Test the get_file function."""

    res = Downloader(directory=tmp_path)
    # Delete file if it already exists, so the test always downloads it
    dl_path = tmp_path / filename
    if dl_path.exists():
        dl_path.unlink()

    if expected is not None:
        with pytest.raises(expected):
            res.get_file(filename)

    else:
        # Download data (or check that it already exists)
        res.get_file(filename)

        # Get the file again, already existing so it doesn't download it again
        res.get_file(filename)


@pytest.mark.flaky(reruns=5)  # in case of intermittent connection to World Wide Web™
def test_get_file_NIST_PSTAR_datafile(tmp_path) -> None:
    """Test the get_file function on a NIST PSTAR datafile."""
    filename = "NIST_PSTAR_aluminum.txt"

    # Delete file if it already exists, so the test always downloads it
    dl_path = tmp_path / filename
    if dl_path.exists():
        dl_path.unlink()

    res = Downloader(directory=tmp_path)
    # Download data (or check that it already exists)
    path = res.get_file(filename)

    arr = np.loadtxt(path, skiprows=7)
    assert np.allclose(arr[0, :], np.array([1e-3, 1.043e2]))


def test_existing_sha_blob_file(tmp_path):
    """
    Test persistence of the blob file between instances
    """
    # Create one object, creating the blob SHA file
    res = Downloader(directory=tmp_path)
    res.blob_dict["test_key"] = "sha_hash"
    res._write_blobfile()

    # Create a second resource object, which should now read the existing
    # blob file
    res = Downloader(directory=tmp_path)

    assert "test_key" in res.blob_dict


def test_multiple_resource_calls(tmp_path):
    """
    Test various file retrieval modes
    """
    # Create a dummy file
    filename = "NIST_PSTAR_aluminum.txt"

    res = Downloader(directory=tmp_path)

    # Download from repository
    res.get_file(filename)

    # Return local copy, since file now exists
    res.get_file(filename)

    # Retrieve a local file that isn't on the remote
    # First create the file
    filename2 = "not_on_the_repo.txt"
    filepath2 = Path(tmp_path, filename2)
    with filepath2.open("w") as f:
        f.write("Not data")
    # Add it to the blob file
    res.blob_dict[filename2] = "sha"
    res._write_blobfile()
    # Now try retrieving it
    res.get_file(filename2)

    # Error is raised when a file isn't local or on the remote
    with pytest.raises(ValueError):
        res.get_file("not_anywhere.txt")

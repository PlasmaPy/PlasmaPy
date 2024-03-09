import numpy as np
import pytest
import os
from pathlib import Path

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


@pytest.mark.parametrize(("filename", "expected"), test_files)
@pytest.mark.flaky(reruns=5)  # in case of intermittent connection to World Wide Web™
def test_get_file(filename, expected, tmp_path) -> None:
    """Test the get_file function."""
    # Delete file if it already exists, so the test always downloads it
    dl_path = tmp_path / filename
    if dl_path.exists():
        dl_path.unlink()

    if expected is not None:
        with pytest.raises(expected):
            downloader.get_file(filename, directory=tmp_path)

    else:
        # Download data (or check that it already exists)
        downloader.get_file(filename, directory=tmp_path)

        # Get the file again, already existing so it doesn't download it again
        downloader.get_file(filename, directory=tmp_path)


@pytest.mark.flaky(reruns=5)  # in case of intermittent connection to World Wide Web™
def test_get_file_NIST_PSTAR_datafile(tmp_path) -> None:
    """Test the get_file function on a NIST PSTAR datafile."""
    filename = "NIST_PSTAR_aluminum.txt"

    # Delete file if it already exists, so the test always downloads it
    dl_path = tmp_path / filename
    if dl_path.exists():
        dl_path.unlink()

    # Download data (or check that it already exists)
    path = downloader.get_file(filename, directory=tmp_path)

    arr = np.loadtxt(path, skiprows=7)
    assert np.allclose(arr[0, :], np.array([1e-3, 1.043e2]))
    
    
    
def test_update_downloads(tmp_path)->None:
    """Test the update_downloads function"""
 
    # Create a file with the same name as a file from the data repository
    # but with different content 
    invalid_data = "Invalid data"
    data_path = Path(tmp_path, "NIST_PSTAR_aluminum.txt")
    with open(data_path, 'w') as f:
        f.write(invalid_data)
        
    # Create another file and directory that shouldn't be touched by the
    # updating function 
    folder_path = Path(tmp_path, 'folder')
    os.makedirs(folder_path, exist_ok=True)

    file_path = Path(tmp_path, 'file1.txt')
    test_data = 'Test data'
    with open(file_path, 'w') as f:
        f.write(test_data)
    
    # Update files in the directory
    files_updated = downloader.update_downloads(directory=tmp_path)
    
    # One file should have been updated 
    assert len(files_updated) == 1
    
    # The first line of the updated file should now be different 
    with open(data_path, 'r') as f:
       first_line = f.readline()
    assert first_line != invalid_data
    
    # The folder and other file should still be there
    assert folder_path.is_dir()
    assert file_path.is_file()
    
    with open(file_path, 'r') as f:
       first_line = f.readline()
    assert first_line == test_data
    

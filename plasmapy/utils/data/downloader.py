"""
Contains functionality for downloading files from a URL. Intended for
downloading files from |PlasmaPy's data repository|.

"""

from pathlib import Path
from urllib.parse import urljoin
from shutil import rmtree

import requests

__all__ = ["get_file", 'update_downloads']

# Note: GitHub links have a problem where the Content-Encoding is set to
# 'gzip' but the file is not actually compressed. This header is just ignored
# by the get_file function.
_BASE_URL = "https://raw.githubusercontent.com/PlasmaPy/PlasmaPy-data/main/"

# TODO: use a config file variable to allow users to set a location
# for the data download folder?

# Default location for downloaded data
_default_downloads_folder = Path(Path.home(), ".plasmapy", "downloads")


def update_downloads(directory:str|None=None)->None:
    """
    Updates all downloaded resource files in the provided directory. 
    
    Parameters
    ----------
    directory : str, optional
        The full path to the desired download location. Defaults to the
        default PlasmaPy data download directory
        :file:`plasmapy/utils/data/downloads/`\ .
        
        
    Returns
    -------
    files_updated : list(str)
        A list of paths for files that were re-downloaded
        
    """
    
    if directory is None:
        directory = _default_downloads_folder
        
    files_updated = []
    for path in directory.glob('**/*'):
        if path.is_file():
            # For each file, try to download again, accepting 404 errors
            # which will occur if the file doesn't exist on the repository
            try:
                path = get_file(path.name, directory=directory, 
                                force_download=True)
                files_updated.append(path)
            except OSError:
                pass
            
    return files_updated
                


def get_file(basename, base_url=_BASE_URL, directory=None, 
             force_download=False):
    r"""
    Download a file from a URL (if the file does not already exist) and
    return the full local path to the file.

    Parameters
    ----------
    basename : str
        Name of the file to be downloaded (extension included).

    base_url : str, optional
        The base URL of the file to be downloaded. Defaults to the root
        directory of |PlasmaPy's data repository|.

    directory : str, optional
        The full path to the desired download location. Defaults to the
        default PlasmaPy data download directory
        :file:`plasmapy/utils/data/downloads/`\ .
        
    force_download : bool, optional
        If True, re-download the file even if it already exists in the
        destination file. The default value is `False`.

    Returns
    -------
    path : str
        The full local path to the downloaded file.
        
        
    Notes
    -----
    Once a file is downloaded, it will not be re-downloaded even if the file
    changes on the remote repository. To update (re-download) all data files,
    use the `~plasmapy.util.downloader.update_downloads` function.
     
    """

    if "." not in str(basename):
        raise ValueError(f"'filename' ({basename}) must include an extension.")

    if directory is None:
        directory = _default_downloads_folder

        # Create the .plasmapy/downloads directory if it does not already
        # exist
        if not directory.is_dir():
            directory.mkdir()

    path = Path(directory, basename)
    
    # If file doesn't exist locally, download it
    if not path.is_file() or force_download:
        url = urljoin(base_url, basename)

        reply = requests.get(url)  # noqa: S113

        if reply.status_code == 404:
            raise OSError(
                "The requested URL returned a 404 code, which "
                "likely indicates that the file does not exist at the "
                "URL provided."
            )
            
        # If force_download is set, remove any existing copy of the file
        # Only do this after the request, in case request raises an 
        # exception 
        if path.is_file():
            path.unlink()

        with path.open(mode="wb") as f:
            f.write(reply.content)

    return path

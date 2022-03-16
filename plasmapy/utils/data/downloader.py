"""
Contains functionality for downloading files from a URL. Intended for
downloading files from the PlasmaPy data repository.

"""

import os
import requests

from urllib.parse import urljoin

__all__ = ["get_file"]

# Note: GitHub links have a problem where the Content-Encoding is set to
# 'gzip' but the file is not actually compressed. This header is just ignored
# by the get_file function.
_BASE_URL = "https://github.com/PlasmaPy/sample-data/raw/main/data/"

# TODO: use a config file variable to allow users to set a location for this folder?
_DOWNLOADS_PATH = os.path.join(os.path.dirname(__file__), "downloads")


def get_file(basename, base_url=_BASE_URL, directory=None):
    r"""
    Downloads a file from a URL (if the file does not already exist) and
    returns the full local path to the file.

    Parameters
    ----------
    basename : str
        Name of the file to be downloaded (extension included).

    base_url : str, optional
        The base URL of the file to be downloaded. Defaults to the main
        directory of the PlasmaPy data repository.

    directory : str, optional
        The full path to the desired download location. Defaults to the
        default PlasmaPy data download directory.

    Returns
    -------
    path : str
        The full local path to the downloaded file.

    """

    if "." not in str(basename):
        raise ValueError(f"'filename' ({basename}) must include an extension.")

    if directory is None:
        directory = _DOWNLOADS_PATH

    path = os.path.join(directory, basename)

    # If file doesn't exist locally, download it
    if not os.path.exists(path):

        url = urljoin(base_url, basename)

        reply = requests.get(url)

        # Missing files on GitHub will resolve to a 404 html page, so we use
        # this as an indicator that the file may not exist.
        if "text/html" in reply.headers["Content-Type"]:
            raise OSError(
                "The requested URL returned an html file, which "
                "likely indicates that the file does not exist at the "
                "URL provided."
            )

        with open(path, "wb") as f:
            f.write(reply.content)

    return path

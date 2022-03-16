"""
Contains functionality for downloading data from the PlasmaPy data repository.

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


def get_file(basename, base_url=_BASE_URL):
    """
    Downloads a file with a given filename from the PlasmaPy/data
    repository.

    Parameters
    ----------
    basename : str
        Name of the file to be downloaded (extension included).

    base_url : str, optional
        The URL of the PlasmaPy data repository.

    Returns
    -------
    path : str
        DESCRIPTION.

    """
    if not "." in str(basename):
        raise ValueError(f"'filename' ({basename}) must include an extension.")

    path = os.path.join(_DOWNLOADS_PATH, basename)

    # If file doesn't exist, download it
    if not os.path.exists(path):

        url = urljoin(base_url, basename)

        # Get the requested content
        r = requests.get(url)

        # Validate that the content type matches one of the content types
        # the module knows how to download.
        #
        # Missing files on GitHub will resolve to a 404 html page, so we use
        # this as an indicator that the file may not exist.
        allowed_types = ["text/plain; charset=utf-8", "image/png"]

        if not r.headers["Content-Type"] in allowed_types:
            raise OSError(
                f"The requested file is not an allowed"
                f"Content-Type: {r.headers['Content-Type']}."
                "This may indicate that the file does not exist at "
                "the URL provided."
            )

        # Write the content to disk
        with open(path, "wb") as f:
            f.write(r.content)

    return path

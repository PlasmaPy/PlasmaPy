"""
Contains functionality for downloading files from a URL. Intended for
downloading files from `PlasmaPy's data repository`_.

"""

import requests

from pathlib import Path
from urllib.parse import urljoin

__all__ = ["get_file"]

# Note: GitHub links have a problem where the Content-Encoding is set to
# 'gzip' but the file is not actually compressed. This header is just ignored
# by the get_file function.
_BASE_URL = "https://raw.githubusercontent.com/PlasmaPy/PlasmaPy-data/main/"

# TODO: use a config file variable to allow users to set a location
# for the data download folder?


def get_file(basename, base_url=_BASE_URL, directory=None):
    r"""
    Download a file from a URL (if the file does not already exist) and
    return the full local path to the file.

    Parameters
    ----------
    basename : str
        Name of the file to be downloaded (extension included).

    base_url : str, optional
        The base URL of the file to be downloaded. Defaults to the root
        directory of `PlasmaPy's data repository`_.

    directory : str, optional
        The full path to the desired download location. Defaults to the
        default PlasmaPy data download directory
        :file:`plasmapy/utils/data/downloads/`\ .

    Returns
    -------
    path : str
        The full local path to the downloaded file.
    """

    if "." not in str(basename):
        raise ValueError(f"'filename' ({basename}) must include an extension.")

    if directory is None:  # coverage: ignore
        directory = Path(Path.home(), ".plasmapy", "downloads")

        # Create the .plasmapy/downloads directory if it does not already
        # exist
        if not directory.is_dir():
            directory.mkdir()

    path = Path(directory, basename)

    # If file doesn't exist locally, download it
    if not path.is_file():
        url = urljoin(base_url, basename)

        reply = requests.get(url)

        if reply.status_code == 404:
            raise OSError(
                "The requested URL returned a 404 code, which "
                "likely indicates that the file does not exist at the "
                "URL provided."
            )

        with path.open(mode="wb") as f:
            f.write(reply.content)

    return path

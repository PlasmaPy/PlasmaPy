"""
Do

"""

import os

from parfive import Downloader
from urllib.parse import urljoin

__all__ = ["get_file"]

# TODO: the normal GitHub raw links seem to have an issue with file compression?
# their Content-Type header seems to be incorrectly set to say they are using gzip, but the files disagree?
# _BASE_URL = "https://github.com/PlasmaPy/sample-data/raw/main/data/"
# _BASE_URL = "https://raw.githubusercontent.com/PlasmaPy/sample-data/main/data/"

# TODO: A more permanent fix for this issue?
# This service fixes the GitHub headers...
_BASE_URL = "https://gitcdn.link/repo/PlasmaPy/data/main/data/"

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
        raise ValueError("'filename' must include an extension.")

    path = os.path.join(_DOWNLOADS_PATH, basename)

    # If file doesn't exist, download it
    if not os.path.exists(path):
        url = urljoin(base_url, basename)
        dl = Downloader(
            overwrite=True, progress=True, headers={"Accept-Encoding": "identity"}
        )
        dl.enqueue_file(url, path=_DOWNLOADS_PATH, filename=basename)
        results = dl.download()

        # If an error is returned, raise an exception
        if len(results.errors) > 0:
            raise OSError(results.errors[0][2])
        # Otherwise, the only element in results is the filepath
        else:
            path = results[0]

    return path


if __name__ == "__main__":
    get_file("ex")

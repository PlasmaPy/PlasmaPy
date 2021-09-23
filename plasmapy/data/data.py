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

# This service fixes the GitHub headers...
_BASE_URL = "https://gitcdn.link/repo/PlasmaPy/sample-data/main/data/"

_DOWNLOADS_PATH = os.path.join(os.path.dirname(__file__), "downloads")


def get_file(filename, base_url=_BASE_URL):
    path = os.path.join(_DOWNLOADS_PATH, filename)

    # If file doesn't exist, download it
    if not os.path.exists(path):
        url = urljoin(base_url, filename)
        dl = Downloader(
            overwrite=True, progress=True, headers={"Accept-Encoding": "identity"}
        )
        dl.enqueue_file(url, path=_DOWNLOADS_PATH, filename=filename)
        files = dl.download()
        path = files[0]
    return path

"""
Contains functionality for downloading files from a URL. Intended for
downloading files from |PlasmaPy's data repository|.

"""

import json
import warnings
from pathlib import Path
from urllib.parse import urljoin

import requests

__all__ = ["Downloader"]


# TODO: use a config file variable to allow users to set a location
# for the data download folder?


class Downloader:
    """
    Accesses the PlasmaPy resource files.

    Retrieves local paths to resource files, and downloads those files from
    |PlasmaPy's data repository| if they cannot be found locally.

    Parameters
    ----------
    directory : Path|None, optional
        DESCRIPTION. The default is None.

    """

    _API_BASE_URL = "https://api.github.com/repos/PlasmaPy/PlasmaPy-data/contents/"

    _blob_file = "RESOURCE_BLOB_SHA.json"

    def __init__(self, directory: Path | None = None):
        if directory is None:
            self._download_directory = Path(
                Path.home(), ".plasmapy", "downloads"
            )  # nocov
        else:
            self._download_directory = Path(directory)

        # Make the directory if it doesn't already exist
        self._download_directory.mkdir(parents=True, exist_ok=True)

        # Path to the SHA blob file
        self._blob_file_path = Path(self._download_directory, self._blob_file)

        # Create the SHA blob file if it doesn't already exist
        if not self._blob_file_path.is_file():
            self.blob_dict = {}
            self._write_blobfile()
        # Otherwise, read the SHA blob file
        else:
            self._read_blobfile()

    def _write_blobfile(self) -> None:
        """
        Write the blob_dict to disk.
        """
        with self._blob_file_path.open("w") as f:
            json.dump(self.blob_dict, fp=f)

    def _read_blobfile(self) -> None:
        """
        Read the blob_dict from disk.
        """
        with self._blob_file_path.open("r") as f:
            self.blob_dict = json.load(f)

    def _http_request(self, url: str) -> requests.Response:
        """
        Issue an HTTP request to the specified URL, handling exceptions.
        """

        reply = requests.get(url)  # noqa: S113

        # Extract the 'message' value if it is there
        # If the file does not exist on the repository, the GitHub API
        # will return `Not Found` in response to this but not raise a 404 error

        if reply.status_code == 404 or b"Not Found" in reply.content:
            raise FileNotFoundError

        return reply

    def _repo_file_info(self, filename: str) -> tuple[str, str]:
        """
        Return file information from github via the API.

        Parameters
        ----------
        filename : str
            DESCRIPTION.

        Returns
        -------
        sha : str
            SHA hash for the file on GitHub.
        dl_url : str
            URL from which the file can be downloaded from GitHub.

        """
        url = urljoin(self._API_BASE_URL, filename)

        try:
            reply = self._http_request(url)
        except FileNotFoundError as err:
            raise FileNotFoundError(f"File {filename} not found at {url}") from err

        # No test coverage for this exception since we can't test it without
        # severing the network connectivity in pytest
        except requests.ConnectionError as err:  # nocov
            raise requests.ConnectionError(  # nocov
                "Unable to connect to data "  # nocov
                f"repository {self._API_BASE_URL}"  # nocov
            ) from err  # nocov

        # Extract the SHA hash and the download URL from the response
        info = reply.json()
        sha = info["sha"]
        dl_url = info["download_url"]

        return sha, dl_url

    def _download_file(self, filepath: Path, dl_url: str) -> None:
        """
        Download a file from a given URL to a specified path.
        """

        # Request the contents of the file from the download URL
        reply = self._http_request(dl_url)

        # Write the contents to file
        with filepath.open(mode="wb") as f:
            f.write(reply.content)

    def get_file(self, filename: str) -> Path:
        """
        Returns a local path to a resource file, downloading it if necessary.

        Parameters
        ----------
        filename : str
            The name of the file in the PlasmaPy-data repository.

        Raises
        ------
        ValueError
            If the file cannot be found locally or online.

        Returns
        -------
        Path : Path
            The local path to the resource file.

        """
        # Update the memory copy of the blob file
        self._read_blobfile()

        filepath = Path(self._download_directory, filename)

        # If local file exists and also exists in blob file, get the
        # file sha
        if filepath.is_file() and filename in self.blob_dict:
            local_sha = self.blob_dict[filename]
        else:
            local_sha = None

        # Get the online SHA
        try:
            online_sha, dl_url = self._repo_file_info(filename)
        # If online file cannot be found, set the sha hash to None
        except (requests.ConnectionError, FileNotFoundError):
            online_sha = None

        # If local sha and online sha are equal, return the local filepath
        if local_sha == online_sha and local_sha is not None:
            return filepath

        # Try downloading from the repository
        elif online_sha is not None:
            # Download the file
            self._download_file(filepath, dl_url)

            # Add SHA to blob dict and update blob file
            self.blob_dict[filename] = online_sha
            self._write_blobfile()

            local_sha = online_sha
            return filepath

        # If online file cannot be reached but local file is present,
        # return local file with warning
        elif online_sha is None and local_sha is not None:
            warnings.warn(
                "Request to PlasmaPy-data repository returned 404: "
                "proceeding with local files only, which may be out "
                "of date."
            )
            return filepath

        # If neither online file or local file can be found, raise an
        # exception
        else:
            raise ValueError(
                "Resource could not be found locally or "
                "retrieved from the PlasmPy-data repository: "
                f"{filename}"
            )

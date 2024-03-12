"""
Contains functionality for downloading files from a URL. Intended for
downloading files from |PlasmaPy's data repository|.

"""

import json
import warnings
from pathlib import Path

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
    directory : `~pathlib.Path`, optional
        The directory into which files will be downloaded. The default
        is /~//.plasmapy//downloads//

    """

    _API_BASE_URL = "https://api.github.com/repos/PlasmaPy/PlasmaPy-data/contents/"

    _local_blob_file = "RESOURCE_BLOB_SHA.json"

    def __init__(self, directory: Path | None = None):
        if directory is None:
            # No test coverage for default directory, since pytest always
            # saves into a temporary directory
            self._download_directory = Path(  # coverage: ignore
                Path.home(),
                ".plasmapy",
                "downloads",
            )
        else:
            self._download_directory = Path(directory)

        # Make the directory if it doesn't already exist
        self._download_directory.mkdir(parents=True, exist_ok=True)

        # Path to the SHA blob file
        self._local_blob_file_path = Path(
            self._download_directory, self._local_blob_file
        )

        # Create the SHA blob file if it doesn't already exist
        if not self._local_blob_file_path.is_file():
            self._local_blob_dict = {}
            self._write_blobfile()
        # Otherwise, read the SHA blob file
        else:
            self._read_blobfile()

        # Retrieve the online repo file information
        try:
            self._repo_blob_dict = self._get_repo_blob_dict()
        except ValueError:
            self._repo_blob_dict = None

    def _write_blobfile(self) -> None:
        """
        Write the _local_blob_dict to disk.
        """
        with self._local_blob_file_path.open("w") as f:
            json.dump(self._local_blob_dict, fp=f)

    def _read_blobfile(self) -> None:
        """
        Read the _local_blob_dict from disk.
        """
        with self._local_blob_file_path.open("r") as f:
            self._local_blob_dict = json.load(f)

    def _http_request(self, url: str) -> requests.Response:
        """
        Issue an HTTP request to the specified URL, handling exceptions.
        """

        try:
            reply = requests.get(url)  # noqa: S113

        # No test coverage for this exception since we can't test it without
        # severing the network connectivity in pytest
        except requests.ConnectionError as err:  # coverage: ignore
            raise requests.ConnectionError(
                "Unable to connect to data " f"repository {self._API_BASE_URL}"
            ) from err

        # Extract the 'message' value if it is there
        # If the file does not exist on the repository, the GitHub API
        # will return `Not Found` in response to this but not raise a 404 error
        if reply.status_code == 404:
            raise ValueError(f"URL returned 404: {url}")

        return reply

    def _get_repo_blob_dict(self) -> dict[dict[str, str]]:
        """
        Download the file information for all the files in the repository
        at once.

        Raises
        ------
        ValueError
            If the URL does not return the expected JSON file with the
            expected keys.

        Returns
        -------
        repo_blob_dict : dict
            Dictionary with filenames as keys. Each item is another entry
            with keys `sha` and `download_url`.
        """

        # Ensure that the GitHub API is not rate limited
        reply = self._http_request("https://api.github.com/rate_limit")
        info = reply.json()
        rate_info = info["resources"]["core"]
        limit = int(rate_info["limit"])
        used = int(rate_info["used"])
        if used >= limit:  # coverage: ignore
            raise ValueError(
                f"Exceeded GitHub API limit ({used}/{limit}), "
                "please try Downloader again later."
            )

        reply = self._http_request(self._API_BASE_URL)

        # Extract the SHA hash and the download URL from the response

        # Extract contents to JSON
        # Not tested, since any URL on the gituhb API that doesn't raise a 404
        # should return a JSON
        try:  # coverage: ignore
            info = reply.json()
        except requests.exceptions.JSONDecodeError as err:
            raise ValueError(
                "URL did not return the expected JSON file: "
                f"{self._API_BASE_URL}. "
                f"Response content: {reply.content}"
            ) from err

        repo_blob_dict = {}

        for item in info:
            try:
                repo_blob_dict[item["name"]] = {
                    "sha": item["sha"],
                    "download_url": item["download_url"],
                }
            # Not tested, since any URL on the gituhb API that doesn't return a 404
            # should be a JSON with these keys
            except KeyError as err:  # coverage: ignore
                raise ValueError(
                    f"URL {self._API_BASE_URL} returned JSON file, "
                    "but missing expected "
                    f"keys 'sha' and 'download_url`. JSON contents: {info}"
                ) from err

            except TypeError as err:
                raise TypeError(f"Unexpected response type {info}") from err

        return repo_blob_dict

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
        Path : `~pathlib.Path`
            The local path to the resource file.

        """
        # Update the memory copy of the blob file
        self._read_blobfile()

        filepath = Path(self._download_directory, filename)

        # If local file exists and also exists in blob file, get the
        # file sha
        if filepath.is_file() and filename in self._local_blob_dict:
            local_sha = self._local_blob_dict[filename]
        else:
            local_sha = None

        # Get the online SHA
        repo_err = None
        online_sha = None
        if self._repo_blob_dict is not None:
            try:
                online_sha = self._repo_blob_dict[filename]["sha"]
            except KeyError as err:
                repo_err = err
                warnings.warn(f"Filename {filename} not found on repository.")
            except ValueError as err:
                repo_err = err

        # If local sha and online sha are equal, return the local filepath
        if local_sha == online_sha and local_sha is not None:
            return filepath

        # Try downloading from the repository
        elif online_sha is not None:
            dl_url = self._repo_blob_dict[filename]["download_url"]
            # Download the file
            self._download_file(filepath, dl_url)

            # Add SHA to blob dict and update blob file
            self._local_blob_dict[filename] = online_sha
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
                f"{filename}. Exception raised: {repo_err}"
            )

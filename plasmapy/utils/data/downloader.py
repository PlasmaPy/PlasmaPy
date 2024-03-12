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
    directory : `~pathlib.Path`, optional
        The directory into which files will be downloaded. The default
        is /~//.plasmapy//downloads//

    validate : bool, optional
        If `True`, verify that local files are up-to-date with the data
        repository, and use the GitHub API to verify download URLs before
        downloading files. If `False`, return any matching local file without
        verification and, if a local file cannot be found, attempt to download
        blindly from the repository. Default is `True`.

    """

    # URL for the PlasmaPy-data repository through the GitHub API
    _API_BASE_URL = "https://api.github.com/repos/PlasmaPy/PlasmaPy-data/contents/"

    # Name for the local file that stores the SHA hashes of downloaded files
    _local_blob_file = "RESOURCE_BLOB_SHA.json"

    # Base URL for RAW files
    _RAW_BASE_URL = "https://raw.githubusercontent.com/PlasmaPy/PlasmaPy-data/main/"

    def __init__(self, directory: Path | None = None, validate: bool = True):
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

        self._validate = validate

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

    def _api_usage(self):
        """
        Return the API call limit and the number currently used from this IP.

        """
        # Ensure that the GitHub API is not rate limited
        reply = self._http_request("https://api.github.com/rate_limit")
        info = reply.json()
        rate_info = info["resources"]["core"]
        limit = int(rate_info["limit"])
        used = int(rate_info["used"])

        return limit, used

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

        # If validation is disabled, do not request any information from
        # the server
        if not self._validate:
            return None

        limit, used = self._api_usage()
        # No tests since we don't want to hit this limit!
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

    def _download_file(self, filename: str, dl_url: str) -> str:
        """
        Download a file from a given URL to a specified path.

        Parameters
        ----------
        filename : str
            Name of the file to download: determines the filepath

        dl_url : str
            URL from which to download


        Returns
        -------
        filepath : str
            Path to the downloaded file

        """

        # Request the contents of the file from the download URL
        reply = self._http_request(dl_url)

        filepath = Path(self._download_directory, filename)

        # Write the contents to file
        with filepath.open(mode="wb") as f:
            f.write(reply.content)

        return filepath

    def _get_file_with_validation(
        self, filename: str, local_sha: str, repo_sha: str
    ) -> str | None:
        """
        Return file logic with validation.
        """

        filepath = Path(self._download_directory, filename)

        # If local sha and online sha are equal, return the local filepath
        if local_sha == repo_sha and local_sha is not None:
            return filepath

        # If the file is found online, try downloading from the repository
        elif repo_sha is not None:
            dl_url = self._repo_blob_dict[filename]["download_url"]
            # Download the file
            filepath = self._download_file(filename, dl_url)

            # Add SHA to blob dict and update blob file
            self._local_blob_dict[filename] = repo_sha
            self._write_blobfile()

            local_sha = repo_sha
            return filepath

        # If online file cannot be reached but local file is present,
        # return local file with warning
        elif repo_sha is None and local_sha is not None:
            warnings.warn(
                "No connection to PlasmaPy-data repository? "
                "Proceeding with local files only, which may be out of date."
            )
            return filepath

        else:
            return None

    def _get_file_without_validation(
        self, filename: str, local_sha: str
    ) -> str | Exception:
        """
        Return file logic without validation.

        Raises
        ------
        ValueError
            If the URL is not valid.
        """
        filepath = Path(self._download_directory, filename)

        # If the file exists locally, return that
        if local_sha is not None:
            return filepath

        # Otherwise try blindly downloading from the base URL
        # Note that downloading directly from the RAW url does not
        # require an API call.
        #
        # Raises a ValueError if the URL is not valid
        try:
            dl_url = urljoin(self._RAW_BASE_URL, filename)
            return self._download_file(filename, dl_url)
        except ValueError as err:
            raise ValueError(f"{filename} was not found at URL {dl_url}") from err

    def _get_local_sha(self, filename: str) -> str | None:
        """
        Get the local file SHA hash.
        """
        filepath = Path(self._download_directory, filename)

        # If local file exists and also exists in blob file, get the
        # file sha
        if filepath.is_file() and filename in self._local_blob_dict:
            return self._local_blob_dict[filename]
        else:
            return None

    def _get_repo_sha(self, filename: str) -> tuple[str | None, Exception | None]:
        """
        Get the online file SHA hash.

        """
        repo_sha = None

        # This is skipped if validate=False, since _repo_blob_dict is then None
        if self._repo_blob_dict is not None:
            try:
                repo_sha = self._repo_blob_dict[filename]["sha"]
            except (KeyError, ValueError):
                warnings.warn(f"Filename {filename} not found on repository.")

        return repo_sha

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

        # Get the SHAs
        local_sha = self._get_local_sha(filename)
        repo_sha = self._get_repo_sha(filename)

        if self._validate:
            filepath = self._get_file_with_validation(filename, local_sha, repo_sha)
        else:
            filepath = self._get_file_without_validation(filename, local_sha)

        if filepath is not None:
            return filepath
        else:
            # If neither online file or local file can be found, raise an
            # exception
            raise ValueError(
                "Resource could not be found locally or "
                "retrieved from the PlasmPy-data repository: "
                f"{filename}. "
            )

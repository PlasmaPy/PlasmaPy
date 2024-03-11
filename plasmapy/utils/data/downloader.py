"""
Contains functionality for downloading files from a URL. Intended for
downloading files from |PlasmaPy's data repository|.

"""

from pathlib import Path
from urllib.parse import urljoin
import json
import requests
import warnings

__all__ = ["Resources"]



# TODO: use a config file variable to allow users to set a location
# for the data download folder?


class Resources():

    _API_BASE_URL = "https://api.github.com/repos/PlasmaPy/PlasmaPy-data/contents/"

    _blob_file = 'RESOURCE_BLOB_SHA.json'
    
    def __init__(self, download_directory:Path|None=None):
        
        if download_directory is None:
            self._download_directory =  Path(Path.home(), ".plasmapy", "downloads")
        else:
            self._download_directory = download_directory
        
        # Create the SHA blob file if it doesn't already exist
        self._blob_file_path = Path(self._download_directory, self._blob_file)
        
        if not self._blob_file_path.is_file():
            self.blob_dict = {}
            self._write_blobfile()
        else:
            self._read_blobfile()
            
    def _write_blobfile(self)->None:
        """
        Write the blob_dict to disk.
        """
        with open(self._blob_file_path, 'w') as f:
            json.dump(self.blob_dict, fp=f)
            
    def _read_blobfile(self)->None:
        """
        Read the blob_dict from disk.
        """
        with open(self._blob_file_path, 'r') as f:
            self.blob_dict = json.load(f)
            
            
    def _http_request(self, url:str)->requests.Response:
        """
        Issue an HTTP request to the specified URL, handling exceptions. 
        """
        reply = requests.get(url) # noqa: S113
        
        if reply.status_code == 404:
            raise OSError(
                "The requested URL returned a 404 code, which "
                "likely indicates that the file does not exist at the "
                "URL provided."
            )
        return reply
        
            
    def _repo_file_info(self, filename:str)->tuple[str,str]:
        """
        Return file information from github via the API

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
        url =  urljoin(self._API_BASE_URL, filename)
        reply = self._http_request(url)
            
        # Extract the SHA hash and the download URL from the response 
        info = reply.json()
        sha = info['sha']
        dl_url = info['download_url']
        
        return sha, dl_url
    
    
    def _download_file(self, filepath:Path, dl_url:str)->None:
        """
        Download a file from a given URL to a specified path. 
        """
        
        # Request the contents of the file from the download URL
        reply = self._http_request(dl_url)
        
        # Write the contents to file
        with filepath.open(mode="wb") as f:
            f.write(reply.content)
        

            
    def get_resource(self, filename:str)->Path:
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
        
        filepath = Path(self._download_directory, filename)
        
        # If local file exists and also exists in blob file, get the
        # file sha
        if filepath.is_file() and filename in self.blob_dict.keys():
            local_sha = self.blob_dict[filename]
        else:
            local_sha = None
            
        # Get the online SHA
        try:
            online_sha, dl_url = self._repo_file_info(filename)
        # If online file cannot be found, set the sha hash to None
        except OSError:
            online_sha = None
            dl_url = None
            
        # If local sha and online sha are equal, return the local filepath
        if local_sha == online_sha and local_sha is not None:
            return filepath
    
        # Try downloading from the repository 
        elif online_sha is not None:
            # Download the file
            self._download_file(filepath, dl_url)
            
            # Add SHA to blob dict and update blob file
            self.blob_dict[filename]=online_sha
            self._write_blobfile()
            
            local_sha = online_sha
            return filepath 
        
        # If online file cannot be reached but local file is present,
        # return local file with warning 
        elif online_sha is None and local_sha is not None:
            warnings.warn("Request to PlasmaPy-data repository returned 404: "
                          "proceding with local files only, which may be out "
                          "of date.")
            return filepath
        
        else:
            raise ValueError("Resource could not be found locally or "
                             "retrieved from the PlasmPy-data repository.")
            
            
        
            
            
            
        
        
            
            
        
            
            
        
    
    
    
    
    


if __name__ == '__main__':

    obj = Resources()
    x = obj.get_resource('NIST_STAR.hdf5')
    print(x)
        
        
        
        
        
        






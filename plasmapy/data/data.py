import os
import pkgutil
from urllib.parse import urljoin
from parfive import Downloader


_BASE_URL = "https://github.com/PlasmaPy/sample-data/raw/main/data/"

#_BASE_URL = "https://raw.githubusercontent.com/PlasmaPy/sample-data/main/data/"


def get_file(filename, base_url = _BASE_URL):
    
    
    
    try:
        path = pkgutil.get_data("plasmapy", os.path.join("data", "downloads", filename))
    except FileNotFoundError:
        # If file doesn't exist, download it
        url = urljoin(base_url, filename)
        print(url)
        dl = Downloader(overwrite=True, progress=True, headers={'Accept-Encoding': 'identity'})
        dl.enqueue_file(url, filename=filename)
        files = dl.download()
        print(files)
        



if __name__ == '__main__':
    file = 'NIST_PSTAR_aluminum.txt'
    
    f = get_file(file)
            
    
    
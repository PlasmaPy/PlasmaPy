""" Useful error messages for optional dependencies that aren't found. """
from typing import Optional

def _template(pkgname: str,
              url: str,
              library: Optional[str] = None) -> ImportError:
    """
    Simple template to prevent text duplication.

    Parameters
    ----------
    pkgname : str
        Package name on pip or conda
    url : str
        Link to install page for given package
    library : str (optional)
        Full package name

    Returns
    -------
    ImportError

    """
    library = pkgname if library is None else library
    template = f"""
    {library} could not be found. Try either

        conda install {pkgname}
        
    on Anaconda environments or
        
        pip install {pkgname}

    in general. In case of trouble refer to
    {url}"""

    return ImportError(template)

h5py_import_error = _template("h5py",
                              "http://docs.h5py.org/en/latest/build.html")

mpl_import_error = _template("matplotlib",
                             "https://matplotlib.org/users/installing.html")

mpmath_import_error = _template("mpmath",
                                "http://mpmath.org/doc/current/setup.html#download-and-installation")

lmfit_import_error = _template("lmfit",
                               "https://lmfit.github.io/lmfit-py/installation.html")


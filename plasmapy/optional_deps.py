""" Useful error messages for optional dependencies that aren't found. """
from typing import Optional


def _optional_import_error_template(
    pkgname: str,
    url: str,
    library: Optional[str] = None,
    conda_channel: Optional[str] = None,
) -> ImportError:
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
    conda_channel: str (optional)
        set to "conda-forge" to add -c conda-forge to the conda install line

    Returns
    -------
    ImportError

    """
    library = pkgname if library is None else library
    conda_channel = f"-c {conda_channel} " if conda_channel is not None else ""

    template = f"""
    {library} could not be found. Try either

        conda install {conda_channel}{pkgname}

    on Anaconda environments or

        pip install {pkgname}

    in general. In case of trouble refer to

        {url}

    (link active as of 2018.10.31 - please report dead links on GitHub!)"""

    return ImportError(template)


h5py_import_error = _optional_import_error_template(
    "h5py", "http://docs.h5py.org/en/latest/build.html"
)

mpl_import_error = _optional_import_error_template(
    "matplotlib", "https://matplotlib.org/users/installing.html"
)

mpmath_import_error = _optional_import_error_template(
    "mpmath", "http://mpmath.org/doc/current/setup.html#download-and-installation"
)

lmfit_import_error = _optional_import_error_template(
    "lmfit", "https://lmfit.github.io/lmfit-py/installation.html"
)

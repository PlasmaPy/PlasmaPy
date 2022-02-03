#!/usr/bin/env python
# https://github.com/pypa/pip/issues/7953#issuecomment-645133255
import numpy as np
import site
import sys

from pathlib import Path
from setuptools import setup
from Cython.Build import cythonize  # has to happen after setup import

site.ENABLE_USER_SITE = "--user" in sys.argv[1:]


def get_cython_files():
    """Find all cython files (*.pyx) in plasmapy."""
    return list(Path().glob("plasmapy/**/*.pyx"))


def build_extensions():
    """Build a list of "extensions" for `cythonize`-ing`"""
    # I'm leaving this here as an example of constructing a list of
    # setuptools.Extension objects that could be passed to cythonize().
    # This form could allow for custom parameters for each extension.
    #
    # extensions = []
    # for cy_file in get_cython_files():
    #     mod_name = ".".join(cy_file.with_suffix("").parts)
    #     extensions.append(Extension(mod_name, [str(cy_file)]))

    extensions = [str(cy_file) for cy_file in get_cython_files()]
    return extensions


if __name__ == "__main__":

    # Get configuration information from all the various subpackages.
    # See the docstring for setup_helpers.update_package_files for more
    # details.
    setup(
        use_scm_version=True,
        zip_safe=False,
        include_dirs=[np.get_include()],
        ext_modules=cythonize(build_extensions(), language_level=3),
    )

    # allow for cleaning of --inplace c and shared object files generated
    # by:
    #     python setup.py build_ext --inplace
    #
    # this get executed with the command
    #       python setup.py clean
    #
    if "clean" in sys.argv:
        from importlib.machinery import EXTENSION_SUFFIXES
        ext_suffixes = EXTENSION_SUFFIXES.copy()
        ext_suffixes.append(".c")

        cy_files = get_cython_files()
        for file in cy_files:
            for ext in ext_suffixes:
                ext_file = file.with_suffix(ext)
                if not ext_file.exists():
                    continue

                print(f"removing '{ext_file}'")

                if "--dry-run" in sys.argv:
                    continue

                ext_file.unlink()  # delete extension file

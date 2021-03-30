.. _plasmapy-install:

*******************
Installing PlasmaPy
*******************

.. note::

   If you would like to contribute to PlasmaPy, please refer to the
   instructions on :ref:`installing PlasmaPy for development
   <install-plasmapy-dev>`.

.. _install-requirements:

Requirements
============

PlasmaPy requires Python version 3.7 or newer.
PlasmaPy requires the following packages for installation:

- `NumPy <https://www.numpy.org/>`_ — 1.18.1 or newer
- `SciPy <https://www.scipy.org/>`_ — 1.2 or newer
- `Astropy <https://www.astropy.org/>`_ — 4.0 or newer
- `pandas <https://pandas.pydata.org/>`_ — 1.0 or newer
- `xarray <http://xarray.pydata.org>`_ — above 0.14
- `tqdm <https://tqdm.github.io/>`_ — 4.56 or newer
- `cached_property <https://pypi.org/project/cached-property/>`_ — 1.5.2 or newer

PlasmaPy also depends on the following packages for optional features:

- `h5py <https://www.h5py.org/>`_ — 2.8 or newer
- `lmfit <https://lmfit.github.io/lmfit-py/>`_ — 1.0.1 or newer
- `matplotlib <https://matplotlib.org/>`_ — 2.0 or newer
- `mpmath <https:htt//mpmath.org/>`_ — 1.0 or newer

.. _install-process:

Installing PlasmaPy
===================


.. _install-pip:

Installation with pip
---------------------

To install the most recent release of PlasmaPy on `PyPI`_
with `pip <https://pip.pypa.io/en/stable/>`_ into an existing Python environment
with both required and optional dependencies, run

.. code:: bash

   pip install plasmapy[all]

PlasmaPy may be installed with the required but not the optional dependencies
with the following command, though this may result in an `ImportError` when
using certain specialized functionality.

.. code:: bash

   pip install plasmapy

.. _install-conda:

.. note::

   In some systems, it may be necessary to use ``pip3`` instead of ``pip``.

Installation with conda
-----------------------

We recommend installing PlasmaPy from a Python environment
created using `Conda`_.  Conda allows us to
create and switch between Python environments that are isolated from
each other and the system installation (in contrast to `this xkcd
<https://xkcd.com/1987/>`_).

After `installing conda <https://conda.io/docs/user-guide/install/>`_,
create a PlasmaPy environment with all required and optional dependencies
by running:

.. code:: bash

    conda create -n env_name python=3.9 plasmapy -c conda-forge

where ``env_name`` is replaced by the name of the environment.
To activate this environment, run:

.. code:: bash

   conda activate env_name

Building and installing from source code
========================================

Obtaining source code
---------------------

Stable release
^^^^^^^^^^^^^^

The source code for the most recent stable release of PlasmaPy can be
downloaded `from PyPI`_ or `from Zenodo`_.

Development version on GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have `git`_ installed on your computer, you may clone
`PlasmaPy's GitHub repository`_ and access source code
from the most recent development version by running:

.. code:: bash

   git clone https://github.com/PlasmaPy/PlasmaPy.git

The repository will be cloned inside a new subdirectory called ``PlasmaPy``.

If you do not have git installed on your computer, then you may download
the most recent source code from `PlasmaPy's GitHub repository`_ by
selecting "Clone or Download", which will give you the option to
download a zip file.

.. note::

   Cloning a repository with HTTPS as above is recommended, but you may
   also `clone a repository using SSH`_ as a more secure alternative.

.. note::

   The :ref:`contributing-to-plasmapy` guide has instructions on how to
   fork a repository and create branches so that you may make pull requests.

Building and installing
-----------------------

In the ``PlasmaPy`` directory, run

.. code:: bash

   pip install -e .[all]

where ``-e`` makes the installation editable and ``[all]`` will ensure that
all optional dependencies are installed.  PlasmaPy could also be installed
by running

.. code:: bash

   python setup.py install

Note, however, that this does not download all the dependencies. Check the
`requirements/requirements.txt` file for the current set.

.. _git: https://git-scm.com/
.. _PlasmaPy's GitHub repository: https://github.com/PlasmaPy/PlasmaPy
.. _Conda: https://conda.io/docs/
.. _PyPI: https://pypi.org/
.. _from PyPI: https://pypi.org/project/plasmapy/
.. _from Zenodo: https://doi.org/10.5281/zenodo.1436011
.. _clone a repository using SSH: https://help.github.com/en/github/using-git/which-remote-url-should-i-use#cloning-with-ssh-urls

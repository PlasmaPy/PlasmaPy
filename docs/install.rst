.. _plasmapy-install:

*******************
Installing PlasmaPy
*******************

.. note::

   If you would like to contribute to PlasmaPy, please refer to the
   instructions on :ref:`installing PlasmaPy for development
   <install-plasmapy-dev>`.

.. contents:: Contents
   :local:

Installing Python
=================

PlasmaPy requires Python_ 3.8 or newer. If you do not have Python_
installed already, please follow these instructions to `download
Python`_ and install it.

.. tip::

   If you have trouble installing `plasmapy` on the most recent version
   of Python_ between October and ∼December (e.g., 3.11 for 2022), try
   installing it on the second most recent version (e.g., 3.10). New
   versions of Python_ are released annually in October, and it can take
   time for the scientific Python ecosystem to catch up.

.. _install-pip:

Installation with pip
=====================

To install the most recent release of `plasmapy` on PyPI_ with pip_ into
an existing Python 3.8+ environment, run:

.. code-block:: bash

   pip install plasmapy

On some systems, it might be necessary to use ``pip3`` instead of
``pip``. For more detailed information, please refer to this tutorial on
`installing packages`_.

.. _install-conda:

Installation with Conda
=======================

`plasmapy` can be installed into Python_ environments created by Conda_.
Conda_ allows us to create and switch between Python_ environments that
are isolated from each other and the system installation.

After `installing Conda`_, `plasmapy` can be installed into an activated
Conda_ environment by running:
.. code-block:: bash

   conda install plasmapy

To install `plasmapy` into another existing Conda_ environment, append
:samp:`-n {env_name}` to the previous command, where :samp:`{env_name}`
is replaced with the name of the environment.

To create a new environment with `plasmapy` installed in it, run:

.. code-block:: bash

    conda create -n env_name -c conda-forge plasmapy

where :samp:`env_name` is replaced by the name of the environment.

To activate this environment, run:

.. code-block:: bash

   conda activate env_name

Installation from source code
=============================

Obtaining source code
---------------------

Official releases
^^^^^^^^^^^^^^^^^

A ZIP_ file containing the source code for official releases of
`plasmapy` can be obtained `from PyPI`_ or `from Zenodo`_.

.. Discuss unzipping here

GitHub repository
^^^^^^^^^^^^^^^^^

If you have git_ installed on your computer, you may clone
`PlasmaPy's GitHub repository`_ and access source code
from the most recent development version by running:

.. code:: bash

   git clone https://github.com/PlasmaPy/PlasmaPy.git

The repository will be cloned inside a new subdirectory called
:file:`PlasmaPy`.

If you do not have git_ installed on your computer, then you may download
the most recent source code from `PlasmaPy's GitHub repository`_ by
going to :guilabel:`Code` and selecting :guilabel:`Download ZIP`.
`Unzipping <https://www.wikihow.com/Unzip-a-File>`__ the file will
create a subdirectory called :file:`PlasmaPy` that contains the source
code.

Building and installing
-----------------------

To install the downloaded version of PlasmaPy, enter the
:file:`PlasmaPy` directory and run:

.. code:: bash

   pip install .

If you expect to make any changes to code within PlasmaPy, instead run:

.. code:: bash

   pip install -e .[developer]

The ``-e`` flag makes the installation editable and ``[developer]``
indicates that all of the dependencies needed for developing PlasmaPy
will be installed.

.. note::

   The :ref:`contributing-to-plasmapy` guide has instructions on how to
   fork a repository and create branches so that you may make
   contributions via pull requests.

.. _Anaconda Navigator: https://www.anaconda.com/products/individual
.. _clone a repository using SSH: https://docs.github.com/en/get-started/getting-started-with-git/about-remote-repositories#cloning-with-ssh-urls
.. _download Python: https://www.python.org/downloads/
.. _from PyPI: https://pypi.org/project/plasmapy
.. _from Zenodo: https://doi.org/10.5281/zenodo.1436011
.. _installing Conda:
.. _installing packages: https://packaging.python.org/en/latest/tutorials/installing-packages/#installing-from-vcs
.. _ZIP: https://en.wikipedia.org/wiki/ZIP_(file_format)

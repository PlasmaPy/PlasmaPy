.. _plasmapy-install:

*******************
Installing PlasmaPy
*******************

.. note::

   - PlasmaPy requires Python_ 3.8 or newer.
   - If you would like to contribute to PlasmaPy, please refer to the
     instructions on :ref:`installing PlasmaPy for development
     <install-plasmapy-dev>`.

.. contents:: Contents
   :local:


.. _install-pip:

Installation with pip
=====================

To install the most recent release of `plasmapy` on PyPI_ with pip_ into
an existing Python 3.8+ environment, run:

.. code:: bash

   pip install plasmapy

.. note::

   On some systems, it may be necessary to use ``pip3`` instead of ``pip``.

.. _install-conda:

Installation with Conda
=======================

Conda_ can also be used to install `plasmapy` into a Python environment.
Conda_ allows us to create and switch between Python environments that
are isolated from each other and the system installation (in contrast to
`this xkcd <https://xkcd.com/1987>`_).

After `installing conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_,
create a PlasmaPy environment with all required and optional dependencies
by running:

.. code:: bash

    conda create -n env_name python plasmapy -c conda-forge

where ``env_name`` is replaced by the name of the environment.
To activate this environment, run:

.. code:: bash

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
   fork a repository and create branches so that you may make pull
   requests.

.. _from PyPI: https://pypi.org/project/plasmapy
.. _from Zenodo: https://doi.org/10.5281/zenodo.1436011
.. _clone a repository using SSH: https://docs.github.com/en/get-started/getting-started-with-git/about-remote-repositories#cloning-with-ssh-urls
.. _ZIP: https://en.wikipedia.org/wiki/ZIP_(file_format)

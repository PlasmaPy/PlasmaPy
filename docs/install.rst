.. _plasmapy-install:

Installation
============
There are multiple options to download the source code for PlasmaPy. The
simplest is to select “Clone or Download” on our `repository on
GitHub`_ page.  This will provide an option to download a zip file plus
information on how to clone the repository. If you have git installed on
your computer and you would like to use HTTPS (which is the default and
easier to set up), then run:

.. code:: bash

   git clone https://github.com/PlasmaPy/PlasmaPy.git

If you have `set up an SSH key`_, then an equivalent and more secure
command is:

.. code:: bash

   git clone git@github.com:PlasmaPy/PlasmaPy.git

The :ref:`contributing-to-plasmapy` guide has instructions on how to
fork a repository so that you may make pull requests.

In the top level directory, run

.. code:: bash

   pip install .

or

.. code:: bash

   python setup.py install

--------------

We are officially `on PyPI`_, so PlasmaPy can be installed via

.. code:: bash

   pip install plasmapy

Though be warned that the version on PyPI is our distribution version
and not suitable for development. If you wish to contribute to the
project, please install from GitHub.

We’re not on conda just yet, but we’re working on it!

.. _repository on GitHub: https://github.com/PlasmaPy/PlasmaPy
.. _set up an SSH key: https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/
.. _on PyPI: https://pypi.org/project/plasmapy/
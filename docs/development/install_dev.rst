.. _install-plasmapy-dev:

***********************************
Installing PlasmaPy for Development
***********************************

Obtaining PlasmaPy source code
==============================

After creating your GitHub account, go to the `PlasmaPy repository on
GitHub <https://github.com/PlasmaPy/plasmapy>`_ and **fork a copy of
PlasmaPy to your account**.

To access Git commands on Windows, try `Git Bash
<https://git-scm.com/downloads>`_.

Next you must **clone your fork to your computer**.  Go to the
directory that will host your PlasmaPy directory, and run one of the
following commands (after changing *your-username* to your username).
If you would like to use HTTPS (which is the default and easier to set
up), then run:

.. code-block:: bash

   git clone https://github.com/your-username/PlasmaPy.git

SSH is a more secure option, but requires you to `set up an SSH key
<https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/>`_
beforehand.  The equivalent SSH command is:

.. code-block:: bash

   git clone git@github.com:your-username/PlasmaPy.git

After cloning, we must tell git where the development version of
PlasmaPy is by running:

.. code-block:: bash

   cd PlasmaPy
   git remote add upstream git://github.com/PlasmaPy/PlasmaPy.git

To check on which remotes exist, run ``git remote -v``.  You should get
something like this:

.. code-block:: bash

   origin	git@github.com:namurphy/PlasmaPy.git (fetch)
   origin	git@github.com:namurphy/PlasmaPy.git (push)
   upstream	git@github.com:PlasmaPy/PlasmaPy.git (fetch)
   upstream	git@github.com:PlasmaPy/PlasmaPy.git (push)

Setting up an environment for development
=========================================

Setup procedures for the two most popular virtual environments, conda
and virtualenv, are listed below.

Conda
-----

To set up a development environment for PlasmaPy, we strongly recommend
the `Anaconda distribution <https://www.anaconda.com/download/>`_.

Activate Anaconda
~~~~~~~~~~~~~~~~~

After installing Anaconda, launch any conda environment. By default,
conda installs a ``root`` environment, which you should be able to
activate via

.. code-block:: bash

   source /home/user/anaconda3/bin/activate root

where ``/home/user/anaconda3/`` can be swapped to wherever your anaconda
installation resides.

On `newer versions of Anaconda <https://conda.io/docs/release-notes
.html#recommended-change-to-enable-conda-in-your-shell>`_ the
recommended activation process has changed to:

.. code-block:: bash

   . /home/user/anaconda3/etc/profile.d/conda.sh
   conda activate

.. note::

   On Windows, the way to do this is via running ``Anaconda Prompt``
   from the Start Menu. ``Git Bash`` may also work if you have added
   Anaconda to ``PATH``.

Create your environment
~~~~~~~~~~~~~~~~~~~~~~~

Having activated Anaconda, enter PlasmaPy's repository root directory
and create an environment with our suggested packages by executing the
following:

.. code-block:: bash

   conda env create -f requirements/environment.yml

You may now enter the environment via

.. code-block:: bash

   source activate plasmapy

.. note::

   On Windows, skip the ``source`` part of the previous command.

In newer Conda versions, the command to run is

.. code-block:: bash

   conda activate plasmapy

Virtualenv
----------

Create a directory for holding the PlasmaPy repository, move into it
and create the virtual environment

.. code-block:: bash

   virtualenv -p python3 .

You may need to make sure that this directory's path doesn't contain
any spaces, otherwise virtualenv may throw an error.

Your virtual environment should now be created. If you run ``ls`` you
will notice that virtualenv has created a number of subdirectories:
``bin/``, ``lib/``, and ``include/``. This is why we're not creating the
virtualenv within the repository itself - so as to not pollute it. To
activate the virtualenv you will run:

.. code-block:: bash

   source ./bin/activate

You should now see that your shell session is prepended with
(plasmapy), like so:

.. code-block:: bash

   (plasmapy) user@name:~/programming/plasmapy$

This indicates that the virtualenv is running. Congratulations!  When
your're done working on PlasmaPy, you can deactivate the virtualenv by
running

.. code-block:: bash

   source deactivate

Now that you have plasmapy on your local computer and you have a
virtual environment, you will want to "install" this development
version of PlasmaPy along with its dependencies. Start by activating
your virtual environment. Next you want install the PlasmaPy
dependencies. One way to do this is to do

.. code-block:: bash

   (plasmapy) user@name:~/programming/plasmapy$ pip install -r requirements/environment.txt

Next, setup the development version of PlasmaPy which you just cloned
by moving into the root directory of the cloned repo and running the
setup.py script there:

.. code-block:: bash

   (plasmapy) user@name:~/programming/plasmapy/PlasmaPy$ pip install -e .


You should now be all set to run development versions of PlasmaPy
modules via ``import PlasmaPy`` in your test scripts!

Running anaconda with virtualenv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are running the Anaconda suite and want to use virtualenv to
setup your virtual environment, you will have to let the system know
where the Python interpreter can be found. On Linux this is done with
(for example, assuming having installed Anaconda into ``~/anaconda3``):

.. code-block:: bash

   export LD_LIBRARY_PATH="$HOME/anaconda3/lib/"

Exporting the library path to the dynamic linker will only last for
the duration of the current shell session.

You will have to add the python library directory to LD_LIBRARY_PATH,
as described in a previous step, prior to activating the virtualenv
for every new shell session.

Installing your own dev version
===============================

To be able to import PlasmaPy from your source version, enter the
repository root and use one of

.. code-block:: bash

   python setup.py develop
   pip install -e .

.. note::

   If you are not working within a virtual environment, this may end in
   a permission error - this can be avoided via also adding the
   ``--user`` flag. But seriously, use a virtual environment and spare
   yourself the trouble.

Either one of these commands will create a soft link to your cloned
repository.  Any changes in Python code you make there will be there
when you ``import plasmapy`` from an interactive session.

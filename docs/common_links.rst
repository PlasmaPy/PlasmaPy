.. These are ReST substitutions and links that can be used throughout the docs
   (and docstrings) because they are added to ``docs/conf.py::rst_epilog``.

.. --------
.. Websites
.. --------



.. ----------------------
.. Nested inline literals
.. ----------------------

.. A workaround for nested inline literals so that the filename will get
   formatted like a file but will be a link. In the text, these get used
   with the syntax for a substitution followed by an underscore to
   indicate that it's for a link: |docs/_static|_

.. For these workarounds, if the replacement is something in single back
   ticks (e.g., `xarray`), then it should also be added to
   nitpick_ignore_regex in docs/conf.py so that it doesn't get counted
   as an error in a nitpicky doc build (e.g., tox -e doc_build_nitpicky).

.. _`astropy.units`: https://docs.astropy.org/en/stable/units/index.html
.. |astropy.units| replace:: `astropy.units`

.. _`CITATION.cff`: https://github.com/PlasmaPy/PlasmaPy/blob/main/CITATION.cff
.. |CITATION.cff| replace:: :file:`CITATION.cff`

.. _git: https://git-scm.com
.. |git| replace:: `git`

.. _h5py: https://www.h5py.org/
.. |h5py| replace:: `h5py`

.. _lmfit: https://lmfit.github.io/lmfit-py/
.. |lmfit| replace:: `lmfit`

.. _mpmath: https://mpmath.org/doc/current/
.. |mpmath| replace:: `mpmath`

.. _nbsphinx: https://nbsphinx.readthedocs.io
.. |nbsphinx| replace:: `nbsphinx`

.. _numba: https://numba.readthedocs.io
.. |numba| replace:: `numba`

.. _pre-commit: https://pre-commit.com
.. |pre-commit| replace:: ``pre-commit``

.. _`.pre-commit-config.yaml`: https://github.com/PlasmaPy/PlasmaPy/blob/main/.pre-commit-config.yaml
.. |.pre-commit-config.yaml| replace:: :file:`.pre-commit-config.yaml`

.. _`pyproject.toml`: https://github.com/PlasmaPy/PlasmaPy/blob/main/pyproject.toml
.. |pyproject.toml| replace:: :file:`pyproject.toml`

.. _`tox.ini`: https://github.com/PlasmaPy/PlasmaPy/blob/main/tox.ini
.. |tox.ini| replace:: :file:`tox.ini`

.. _xarray: https://docs.xarray.dev
.. |xarray| replace:: `xarray`

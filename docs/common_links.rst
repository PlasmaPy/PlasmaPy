.. These are ReST substitutions and links that can be used throughout the docs
   (and docstrings) because they are added to ``docs/conf.py::rst_epilog``.

.. --------
.. Websites
.. --------

.. _Astropy: https://docs.astropy.org
.. _BibTeX format: https://www.bibtex.com/g/bibtex-format
.. _black: https://black.readthedocs.io
.. _Contributor Covenant: https://www.contributor-covenant.org
.. _Citation File Format: https://citation-file-format.github.io/
.. _create an issue: https://github.com/PlasmaPy/PlasmaPy/issues/new/choose
.. _CSS: https://www.w3schools.com:443/css
.. _Cython: https://cython.org/
.. _DOI: https://www.doi.org
.. _editable installation: https://pip.pypa.io/en/stable/topics/local-project-installs/#editable-installs
.. _equivalencies: https://docs.astropy.org/en/stable/units/equivalencies.html
.. _flake8: https://flake8.pycqa.org/en/latest
.. _GitHub Actions: https://docs.github.com/en/actions
.. _GitHub Discussions page: https://github.com/PlasmaPy/PlasmaPy/discussions
.. _GitHub Flavored Markdown: https://github.github.com/gfm
.. _GitHub: https://github.com
.. _Gitter bridge: https://gitter.im/PlasmaPy/Lobby
.. _Graphviz: https://graphviz.org
.. _hypothesis: https://hypothesis.readthedocs.io
.. _intersphinx: https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html
.. _isort: https://pycqa.github.io/isort
.. _Jinja: https://jinja.palletsprojects.com
.. _Jupyter: https://jupyter.org
.. _LaTeX: https://www.latex-project.org
.. _mailing list: https://groups.google.com/forum/#!forum/plasmapy
.. _make: https://www.gnu.org/software/make
.. _Markdown: https://www.markdownguide.org
.. _MathJax: https://www.mathjax.org
.. _matplotlib: https://matplotlib.org
.. _Matrix chat room: https://app.element.io/#/room/#plasmapy:openastronomy.org
.. _nbqa: https://nbqa.readthedocs.io
.. _numpydoc: https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard
.. _NumPy: https://numpy.org
.. _office hours: https://www.plasmapy.org/meetings/office_hours/
.. _OpenPMD: https://www.openpmd.org/
.. _pandas: https://pandas.pydata.org
.. _pip: https://pip.pypa.io
.. _Plasma Hack Week: https://hack.plasmapy.org
.. _PlasmaPy: https://www.plasmapy.org
.. _PlasmaPy meetings: https://www.plasmapy.org/meetings
.. _PlasmaPy's documentation: https://docs.plasmapy.org/en/stable
.. _PlasmaPy's GitHub repository: https://github.com/PlasmaPy/plasmapy
.. _PlasmaPy's data repository: https://github.com/PlasmaPy/PlasmaPy-data
.. _PlasmaPy's Matrix chat room: https://app.element.io/#/room/#plasmapy:openastronomy.org
.. _`pre-commit.ci`: https://pre-commit.ci
.. _pydocstyle: https://www.pydocstyle.org/en/stable
.. _pygments: https://pygments.org
.. _PyPI: https://pypi.org
.. _pytest: https://docs.pytest.org
.. _Python: https://www.python.org
.. _Python's documentation: https://docs.python.org/3
.. _Read the Docs: https://readthedocs.org
.. _reST: https://docutils.sourceforge.io/rst.html
.. _reStructuredText (reST): https://docutils.sourceforge.io/rst.html
.. _ruff: https://beta.ruff.rs/docs
.. _SciPy: https://scipy.org
.. _sphinx_automodapi: https://sphinx-automodapi.readthedocs.io
.. _sphinx-build: https://www.sphinx-doc.org/en/master/man/sphinx-build.html
.. _Sphinx: https://www.sphinx-doc.org
.. _suggestion box: https://docs.google.com/forms/d/e/1FAIpQLSdT3O5iHZrLJRuavFyzoR23PGy0Prfzx2SQOcwJGWtvHyT2lw/viewform?usp=sf_link
.. _towncrier: https://github.com/twisted/towncrier
.. _tox: https://tox.wiki/en/latest
.. _virtualenv: https://pypi.org/project/virtualenv
.. _weekly tests: https://github.com/PlasmaPy/PlasmaPy/actions/workflows/weekly.yml
.. _Wikipedia: https://www.wikipedia.org
.. _Zenodo: https://zenodo.org

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

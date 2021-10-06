.. These are ReST substitutions and links that can be used throughout the docs
.. (and docstrings) because they are added to ``docs/conf.py::rst_epilog``.

.. ------------------
.. plasmapy.formulary
.. ------------------

.. |ClassicalTransport| replace:: :class:`~plasmapy.formulary.braginskii.ClassicalTransport`

.. ------------------
.. plasmapy.particles
.. ------------------

.. |CustomParticle| replace:: :class:`~plasmapy.particles.particle_class.CustomParticle`
.. |DimensionlessParticle| replace:: :class:`~plasmapy.particles.particle_class.DimensionlessParticle`
.. |IonicLevel| replace:: :class:`~plasmapy.particles.ionization_state.IonicLevel`
.. |IonizationState| replace:: :class:`~plasmapy.particles.ionization_state.IonizationState`
.. |IonizationStateCollection| replace:: :class:`~plasmapy.particles.ionization_state_collection.IonizationStateCollection`
.. |Particle| replace:: :class:`~plasmapy.particles.particle_class.Particle`
.. |particle_input| replace:: :func:`~plasmapy.particles.decorators.particle_input`
.. |ParticleLike| replace:: :obj:`~plasmapy.particles.particle_class.ParticleLike`
.. |ParticleList| replace:: :class:`~plasmapy.particles.particle_collections.ParticleList`

.. |ChargeError| replace:: :class:`~plasmapy.particles.exceptions.ChargeError`
.. |InvalidElementError| replace:: :class:`~plasmapy.particles.exceptions.InvalidElementError`
.. |InvalidIonError| replace:: :class:`~plasmapy.particles.exceptions.InvalidIonError`
.. |InvalidIsotopeError| replace:: :class:`~plasmapy.particles.exceptions.InvalidIsotopeError`
.. |InvalidParticleError| replace:: :class:`~plasmapy.particles.exceptions.InvalidParticleError`
.. |MissingParticleDataError| replace:: :class:`~plasmapy.particles.exceptions.MissingParticleDataError`
.. |MissingParticleDataWarning| replace:: :class:`~plasmapy.particles.exceptions.MissingParticleDataWarning`
.. |ParticleError| replace:: :class:`~plasmapy.particles.exceptions.ParticleError`
.. |ParticleWarning| replace:: :class:`~plasmapy.particles.exceptions.ParticleWarning`
.. |UnexpectedParticleError| replace:: :class:`~plasmapy.particles.exceptions.UnexpectedParticleError`

.. ------------------
.. NumPy replacements
.. ------------------

.. |ndarray| replace:: :class:`numpy.ndarray`

.. --------------------
.. Astropy replacements
.. --------------------

.. |Quantity| replace:: :class:`~astropy.units.Quantity`
.. |Time| replace:: :class:`~astropy.time.Time`
.. |TimeDelta| replace:: :class:`~astropy.time.TimeDelta`
.. |Unit| replace:: :class:`~astropy.units.UnitBase`

.. --------
.. Websites
.. --------

.. _Astropy docs: https://docs.astropy.org
.. _Astropy: https://www.astropy.org
.. _BibTeX: http://www.bibtex.org/
.. _BibTeX format: https://www.bibtex.com/g/bibtex-format
.. _Conda: https://conda.io/en/latest
.. _CSS: https://en.wikipedia.org/wiki/CSS
.. _docstring: https://en.wikipedia.org/wiki/Docstring
.. _DOI: https://www.doi.org/
.. _GitHub Actions: https://docs.github.com/en/actions
.. _GitHub Discussions page: https://github.com/PlasmaPy/PlasmaPy/discussions
.. _GitHub: https://github.com
.. _Gitter bridge: https://gitter.im/PlasmaPy/Lobby
.. _Jinja: https://jinja.palletsprojects.com
.. _Jupyter: https://jupyter.org
.. _LaTeX: https://www.latex-project.org
.. _Markdown: https://www.markdownguide.org
.. _MathJax: https://www.mathjax.org
.. _Matrix chat room: https://app.element.io/#/room/#plasmapy:openastronomy.org
.. _NumPy: https://numpy.org
.. _pandas: https://pandas.pydata.org
.. _Plasma Hack Week: https://hack.plasmapy.org
.. _PlasmaPy meetings: https://www.plasmapy.org/meetings
.. _PlasmaPy's GitHub repository: https://github.com/PlasmaPy/plasmapy
.. _PlasmaPy's documentation: https://docs.plasmapy.org/en/stable
.. _PlasmaPy: https://www.plasmapy.org
.. _PlasmaPy's Matrix chat room: https://app.element.io/#/room/#plasmapy:openastronomy.org
.. _PyPI: https://pypi.org
.. _Python's documentation: https://docs.python.org/3
.. _Python: https://www.python.org
.. _Read the Docs: https://readthedocs.org
.. _SciPy: https://www.scipy.org
.. _Sphinx: https://www.sphinx-doc.org
.. _Zenodo: https://zenodo.org
.. _black: https://black.readthedocs.io
.. _docs/common_links.rst: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/common_links.rst
.. _git: https://git-scm.com
.. _GitHub Flavored Markdown: https://github.github.com/gfm
.. _Graphviz: https://graphviz.org
.. _intersphinx: https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html
.. _isort: https://pycqa.github.io/isort
.. _mailing list: https://groups.google.com/forum/#!forum/plasmapy
.. _make: https://www.gnu.org/software/make
.. _matplotlib: https://matplotlib.org
.. _numpydoc: https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard
.. _persistent identifier: https://en.wikipedia.org/wiki/Persistent_identifier
.. _pip: https://pip.pypa.io
.. _reST: https://docutils.sourceforge.io/rst.html
.. _reStructuredText (reST): https://docutils.sourceforge.io/rst.html
.. _sphinx_automodapi: https://sphinx-automodapi.readthedocs.io
.. _sphinx-build: https://www.sphinx-doc.org/en/master/man/sphinx-build.html
.. _suggestion box: https://docs.google.com/forms/d/e/1FAIpQLSdT3O5iHZrLJRuavFyzoR23PGy0Prfzx2SQOcwJGWtvHyT2lw/viewform?usp=sf_link
.. _towncrier: https://towncrier.readthedocs.io/en/actual-freaking-docs
.. _tox: https://tox.readthedocs.io
.. _virtualenv: https://pypi.org/project/virtualenv
.. _Wikipedia: https://www.wikipedia.org

.. A workaround for nested inline literals so that the filename will get
   formatted like a file but will be a link. In the text, these get used
   with the syntax for a substitution followed by an underscore to
   indicate that it's for a link: |docs/_static|_

.. _`docs/_static`: https://github.com/PlasmaPy/PlasmaPy/tree/main/docs/_static
.. |docs/_static| replace:: :file:`docs/_static`

.. _`docs/_static/sphinx_rtd_overrides.css`: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/_static/rtd_theme_overrides.css
.. |docs/_static/sphinx_rtd_overrides.css| replace:: :file:`docs/_static/sphinx_rtd_overrides.css`

.. _`docs/api_static`: https://github.com/PlasmaPy/PlasmaPy/tree/main/docs/api_static
.. |docs/api_static| replace:: :file:`docs/api_static`

.. _`docs/conf.py`: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/conf.py
.. |docs/conf.py| replace:: :file:`docs/conf.py`

.. _`docs/glossary.rst`: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/glossary.rst
.. |docs/glossary.rst| replace:: :file:`docs/glossary.rst`

.. _`docs/common_links.rst`: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/common_links.rst
.. |docs/common_links.rst| replace:: :file:`docs/common_links.rst`

.. _`docs/bibliography.bib`: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/bibliography.bib
.. |docs/bibliography.bib| replace:: :file:`docs/bibliography.bib`

.. _`IPython.sphinxext.ipython_console_highlighting`: https://ipython.readthedocs.io/en/stable/sphinxext.html?highlight=IPython.sphinxext.ipython_console_highlighting#ipython-sphinx-directive-module
.. |IPython.sphinxext.ipython_console_highlighting| replace:: ``IPython.sphinxext.ipython_console_highlighting``

.. _nbsphinx: https://nbsphinx.readthedocs.io
.. |nbsphinx| replace:: `nbsphinx`

.. _`setup.cfg`: https://github.com/PlasmaPy/PlasmaPy/blob/main/setup.cfg
.. |setup.cfg| replace:: :file:`setup.cfg`

.. _`sphinxcontrib-bibtex`: https://sphinxcontrib-bibtex.readthedocs.io
.. |sphinxcontrib-bibtex| replace:: `sphinxcontrib-bibtex`

.. _`sphinx_copybutton`: https://sphinx-copybutton.readthedocs.io
.. |sphinx_copybutton| replace:: `sphinx_copybutton`

.. _`sphinx_gallery.load_style`: https://sphinx-gallery.github.io/stable/advanced.html?highlight=load_style#using-only-sphinx-gallery-styles
.. |sphinx_gallery.load_style| replace:: `sphinx_gallery.load_style`

.. _`sphinx_changelog`: https://sphinx-changelog.readthedocs.io
.. |sphinx_changelog| replace:: `sphinx_changelog`

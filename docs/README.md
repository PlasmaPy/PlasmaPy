# PlasmaPy Documentation

> [!TIP]
> To learn more about PlasmaPy's documentation, please check out the
> [**documentation guide**].

PlasmaPy's documentation is written in [reStructuredText], built using
[Sphinx] via a [Nox] session, and hosted by [Read the Docs].

The "stable" documentation corresponds to PlasmaPy's most recent
release, and can be found at https://docs.plasmapy.org/en/stable. The
"latest" documentation corresponds to the current state of the `main`
branch on PlasmaPy's GitHub repository, and is available at
https://docs.plasmapy.org/en/latest.

The [`docs/`] directory contains the source files for PlasmaPy's
narrative documentation. The configuration file is [`docs/conf.py`].

## Building documentation

> [!TIP]
> When making a pull request, the documentation can be previewed by
> clicking on *Details* next to **docs/readthedocs.org:plasmapy** in the
> list of checks.

Building PlasmaPy's documentation requires using the most recent version
of Python supported by PlasmaPy. Prior to building documentation
locally, please install [Nox] and its dependencies with:

```shell
python -m pip install nox uv
```

> [!NOTE]
> It may also be necessary to [install graphviz] and [install pandoc].

The documentation can be built by going to the top-level directory of
your clone of PlasmaPy and running:

```shell
nox -s docs
```

The documentation preview will be built in `docs/build/html/` in your
local clone of PlasmaPy.

[**documentation guide**]: https://docs.plasmapy.org/en/latest/contributing/doc_guide.html
[install graphviz]: https://graphviz.org/download
[install pandoc]: https://pandoc.org/installing.html
[nox]: https://nox.thea.codes
[read the docs]: https://about.readthedocs.com
[restructuredtext]: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#rst-primer
[sphinx]: https://www.sphinx-doc.org
[`docs/conf.py`]: conf.py
[`docs/`]: .

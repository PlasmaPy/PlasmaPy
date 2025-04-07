# Type stubs

This directory contains [type stub files] (ending in `.pyi`) that
provide type information for external packages which have not
implemented [type hint annotations] in their source code. Including
external type information here enables tools like [mypy] to perform
static type checking on code in PlasmaPy that make use of functionality
from these packages.

[mypy]: https://mypy.readthedocs.io
[type hint annotations]: https://docs.python.org/3/glossary.html#term-type-hint
[type stub files]: https://mypy.readthedocs.io/en/stable/stubs.html

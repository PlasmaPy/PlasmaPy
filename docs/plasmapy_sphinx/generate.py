"""
This module contains functionality for auto-generating the stub files related to
the :rst:dir:`automodapi` and :rst:dir:`automodsumm` directives.
"""
__all__ = ["AutomodsummEntry"]

from sphinx.ext.autosummary.generate import (
    AutosummaryEntry,
    AutosummaryRenderer,
    generate_autosummary_content,
)


class AutomodsummEntry(AutosummaryEntry):
    """
    A typed version of `~collections.namedtuple` representing an stub file
    entry for :rst:dir:`automodsumm`.
    """
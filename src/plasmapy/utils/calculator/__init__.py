"""
Subpackage formerly containing a prototype plasma calculator.

.. important::

   The plasma calculator has been removed from PlasmaPy
   so that it can be included in the |plasmapy-calculator| standalone
   package.
"""

from typing import NoReturn  # coverage: ignore


def main() -> NoReturn:  # coverage: ignore
    """Plasma calculator."""
    raise RuntimeError(
        "The plasma calculator has been extracted from PlasmaPy into a "
        "standalone package. For more details, see: "
        "https://github.com/PlasmaPy/plasmapy-calculator"
    )

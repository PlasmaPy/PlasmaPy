import pytest

from plasmapy.utils.decorators.deprecation import deprecated
from plasmapy.utils.exceptions import PlasmaPyDeprecationWarning


def test_deprecated() -> None:
    @deprecated(since="0.7.0")
    def deprecated_function() -> None:
        pass

    with pytest.warns(PlasmaPyDeprecationWarning):
        deprecated_function()

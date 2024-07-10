from datetime import date

import pytest

from plasmapy import __version__


@pytest.mark.skipif("dev" not in __version__, reason="test is irrelevant for ")
def test_dev_version():
    today = date.today()  # noqa: DTZ011
    year = today.year
    month = today.month
    assert __version__.startswith(f"{year}.{month}.")

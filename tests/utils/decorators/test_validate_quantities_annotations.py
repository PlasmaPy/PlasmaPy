"""Test that @validate_quantities works with postponed evaluation of annotations."""

from __future__ import annotations

import astropy.units as u

from plasmapy.utils.decorators.validators import validate_quantities


@validate_quantities  # type: ignore[misc]
def annotated_function(mass: u.Quantity[u.g]) -> u.Quantity[u.kg]:
    return mass


def test_validate_quantities_postponed_annotations() -> None:
    result = annotated_function(1 * u.g)
    assert result.unit == u.kg

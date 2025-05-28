"""Helper functionality for `astropy.units`."""

from __future__ import annotations

__all__: list[str] = []

from numbers import Number
from typing import TYPE_CHECKING

import astropy.units as u

if TYPE_CHECKING:
    from collections.abc import Iterable


def _get_physical_type_dict(
    iterable: Iterable[object],
    *,
    only_quantities: bool | None = False,
    numbers_become_quantities: bool | None = False,
    strict: bool | None = False,
    allowed_physical_types: set[str | u.PhysicalType] | None = None,
) -> dict[u.PhysicalType, u.Quantity]:
    """
    Return a `dict` that contains `~astropy.units.PhysicalType` objects
    as keys and the corresponding objects in ``iterable`` as values.

    Objects in ``iterable`` that do not correspond to a |PhysicalType|
    are skipped.

    Parameters
    ----------
    iterable : iterable
        A iterable that is expected to contain objects with physical
        types.

    only_quantities : `bool`, |keyword-only|, default: `False:
        If `True`, only `~astropy.units.Quantity` instances in
        ``iterable`` will be passed into the resulting `dict`. If
        `False`, then any unit, |PhysicalType|, or object that can be
        converted to a |Quantity| or that has a physical type will be
        included in the `dict`.

    numbers_become_quantities : `bool`, |keyword-only|, default: `False`
        If `True`, numbers will be converted into dimensionless
        |Quantity| instances. If `False`, numbers will be skipped.

    strict : `bool`, |keyword-only|, default: False
        If `True`, raise a `TypeError` if ``iterable`` provides an
        object that does not have a physical type.

    allowed_physical_types : `set` of `~astropy.units.PhysicalType`
        If provided, then if any objects provided by ``iterable`` do not
        have a physical type in ``allowed_physical_types``, then a
        `ValueError` will be raised.

    Returns
    -------
    physical_types : `dict`
        A mapping from |PhysicalType| instances to the corresponding
        objects in ``iterable``.

    Examples
    --------
    >>> import astropy.units as u
    >>> from plasmapy.utils._units_helpers import _get_physical_type_dict
    >>> quantities = [1 * u.m, 2 * u.kg]
    >>> _get_physical_type_dict(quantities)
    {PhysicalType('length'): <Quantity 1. m>, PhysicalType('mass'): <Quantity 2. kg>}
    """
    physical_types = {}

    for obj in iterable:
        if isinstance(obj, Number) and numbers_become_quantities:
            obj = u.Quantity(obj, u.dimensionless_unscaled)  # noqa: PLW2901

        if only_quantities and not isinstance(obj, u.Quantity):
            if strict:
                raise TypeError(f"{obj} is not a Quantity, but should be.")
            continue

        try:
            physical_type = u.get_physical_type(obj)
        except (TypeError, ValueError) as exc:
            if strict:
                raise TypeError(f"{obj} does not have a physical type.") from exc
        else:
            if allowed_physical_types and physical_type not in allowed_physical_types:
                raise ValueError(
                    f"{obj} has a physical type of {physical_type}, but "
                    "only the following physical types are allowed: "
                    f"{allowed_physical_types}."
                )
            if physical_type in physical_types:
                raise ValueError(f"Duplicate physical type: {physical_type}")
            physical_types[physical_type] = obj

    return physical_types

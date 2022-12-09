"""Helper functionality for `astropy.units`."""

from __future__ import annotations

__all__ = []

import astropy.units as u

from numbers import Number
from typing import Dict, Iterable


def _get_physical_type_dict(
    iterable: Iterable,
    *,
    only_quantities=False,
    numbers_become_quantities=False,
) -> Dict[u.PhysicalType, u.Quantity]:
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

    only_quantities : `bool`, |keyword-only|, optional
        If `True`, only `~astropy.units.Quantity` instances in
        ``iterable`` will be passed into the resulting `dict`. If
        `False`, then any unit, |PhysicalType|, or object that can be
        converted to a |Quantity| or that has a physical type will be
        included in the `dict`. Defaults to `False`.

    numbers_become_quantities : `bool`, |keyword-only|, optional
        If `True`, `~numbers.Number` objects will be converted into
        dimensionless |Quantity| instances. If `False`,
        `~numbers.Number` objects will be skipped. Defaults to `False`.

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
            obj = u.Quantity(obj, u.dimensionless_unscaled)

        if only_quantities and not isinstance(obj, u.Quantity):
            continue

        try:
            physical_type = u.get_physical_type(obj)
        except (TypeError, ValueError):
            pass
        else:
            if physical_type in physical_types:
                raise ValueError(f"Duplicate physical type: {physical_type}")
            physical_types[physical_type] = obj

    return physical_types

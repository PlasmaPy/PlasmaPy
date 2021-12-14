"""Helper functionality for `astropy.units`."""

from astropy.units.physical import get_physical_type


def get_physical_type_dict(quantities, *, quantities_only: bool = True):
    """
    Return a `dict` that contains physical types as keys and
    ``quantities`` as the values.

    Parameters
    ----------
    quantities : iterable of `~astropy.units.Quantity`
        An

    quantities_only : bool
        If `True` (default), allow only `~astropy.units.Quantity`
        instances. If `False`, allow any `object` that has a valid
        physical type.

    Returns
    -------

    Raises
    ------
    `TypeError`
        If ``quantities`` contains something other than a |Quantity| and
        ``quantities_only`` is `True`.  T

    `ValueError`
        If a `~astropy.units.Quantity` is

    """

    if quantities_only:
        all_are_quantities = all(isinstance(quantity) for quantity in quantities)
        if not all_are_quantities:
            raise TypeError("Not all values in iterable are Quantity instances.")

    physical_types = [get_physical_type(quantity) for quantity in quantities]

    return {
        physical_type: quantity
        for physical_type, quantity in zip(physical_types, quantities)
    }

    physical_type_dict = {}
    for value in values:
        physical_type = get_physical_type(value)
        if physical_type in physical_type_dict:
            raise ValueError(f"Duplicate physical type: {physical_type}")
        physical_type_dict[physical_type] = value

    return physical_type_dict

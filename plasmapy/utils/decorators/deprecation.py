"""Decorators to mark objects that are deprecated."""

__all__ = ["deprecated"]

import inspect

from astropy.utils.decorators import deprecated as astropy_deprecated

from plasmapy.utils.exceptions import PlasmaPyDeprecationWarning


def deprecated(*args, warning_type=PlasmaPyDeprecationWarning, **kwargs):
    """
    A wrapper of `astropy.utils.decorators.deprecated` that by default assumes
    a warning type of `~plasmapy.utils.exceptions.PlasmaPyDeprecationWarning`.
    """
    return astropy_deprecated(*args, warning_type=warning_type, **kwargs)


# override deprecated's signature
asig = inspect.signature(astropy_deprecated)
new_default = asig.parameters["warning_type"].replace(
    default=PlasmaPyDeprecationWarning
)
new_sig = dict(asig.parameters)
new_sig["warning_type"] = new_default
new_sig = inspect.Signature(
    parameters=tuple(new_sig.values()),
    return_annotation=asig.return_annotation,
)
deprecated.__signature__ = new_sig

"""
Decorators to mark objects as deprecated.

.. deprecated::

   This module has been deprecated in favor of ``@warnings.deprecated``,
   and will be removed in a forthcoming release of PlasmaPy.
"""

__all__ = ["deprecated"]

import inspect

from astropy.utils.decorators import deprecated as astropy_deprecated

from plasmapy.utils.exceptions import PlasmaPyDeprecationWarning


def deprecated(*args, warning_type=PlasmaPyDeprecationWarning, **kwargs):
    """
    A wrapper of `astropy.utils.decorators.deprecated` that by default assumes
    a warning type of `~plasmapy.utils.exceptions.PlasmaPyDeprecationWarning`.

    .. deprecated::

       The decorator has been deprecated in favor of ``@warnings.decorated``,
       and will be removed in a forthcoming release of PlasmaPy.
    """
    return astropy_deprecated(*args, warning_type=warning_type, **kwargs)


# override deprecated's signature
asig = inspect.signature(astropy_deprecated)
new_default = asig.parameters["warning_type"].replace(
    default=PlasmaPyDeprecationWarning,
)
new_sig = dict(asig.parameters)
new_sig["warning_type"] = new_default
new_sig = inspect.Signature(
    parameters=tuple(new_sig.values()),
    return_annotation=asig.return_annotation,
)
deprecated.__signature__ = new_sig

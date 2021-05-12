__all__ = ["deprecated"]

import astropy

from functools import partial

from plasmapy.utils.exceptions import PlasmaPyDeprecationWarning

deprecated = partial(
    astropy.utils.decorators.deprecated,
    warning_type=PlasmaPyDeprecationWarning,
)
"""
A wrapper of `astropy.utils.decorators.deprecated` that by default assumes
a warning type of `~plasmapy.utils.exceptions.PlasmaPyDeprecationWarning`.
"""

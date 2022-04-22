"""
Decorator to convert units of functions in /physics methods
"""
__all__ = ["angular_freq_to_hz"]

import astropy.units as u
import functools
import inspect

from plasmapy.utils.decorators.helpers import preserve_signature


def angular_freq_to_hz(fn):
    """
    A decorator that adds to a function the ability to convert the function's return from
    angular frequency (rad/s) to frequency (Hz).

    A kwarg ``to_hz`` is added to the function's signature, with a default value of `False`.
    The keyword is also added to the function's docstring under the **"Other Parameters"**
    heading.

    Parameters
    ----------
    fn : function
        The function to be decorated

    Raises
    ------
    ValueError
        If ``fn`` has already defined a kwarg ``to_hz``

    Returns
    -------
    callable
        The decorated function

    Notes
    -----
    * If `~plasmapy.utils.decorators.converter.angular_freq_to_hz` is used with decorator
      :func:`~plasmapy.utils.decorators.validators.validate_quantities`, then
      `angular_freq_to_hz` should be used inside
      :func:`~plasmapy.utils.decorators.validators.validate_quantities` but special
      consideration is needed for setup.  The following is an example of an
      appropriate setup::

        import astropy.units as u
        from plasmapy.utils.decorators.converter import angular_freq_to_hz
        from plasmapy.utils.decorators.validators import validate_quantities

        @validate_quantities(validations_on_return={'units': [u.rad / u.s, u.Hz]})
        @angular_freq_to_hz
        def foo(x: u.rad / u.s) -> u.rad / u.s
            return x

      Adding ``u.Hz`` to the allowed units allows the converted quantity to pass
      the validations.

    Examples
    --------

        >>> import astropy.units as u
        >>> from plasmapy.utils.decorators.converter import angular_freq_to_hz
        >>>
        >>> @angular_freq_to_hz
        ... def foo(x):
        ...     return x
        >>>
        >>> foo(5 * u.rad / u.s, to_hz=True)
        <Quantity 0.79577472 Hz>
        >>>
        >>> foo(-1 * u.rad / u.s, to_hz=True)
        <Quantity -0.15915494 Hz>

    Decoration also works with methods

        >>> class Foo:
        ...     def __init__(self, x):
        ...         self.x = x
        ...
        ...     @angular_freq_to_hz
        ...     def bar(self):
        ...         return self.x
        >>>
        >>> foo = Foo(0.5 * u.rad / u.s)
        >>> foo.bar(to_hz=True)
        <Quantity 0.07957747 Hz>

    """
    # raise exception if fn uses the 'to_hz' kwarg
    sig = inspect.signature(fn)
    if "to_hz" in sig.parameters:
        raise ValueError(
            f"Wrapped function '{fn.__name__}' can not use keyword 'to_hz'."
            f" Keyword reserved for decorator functionality."
        )

    # make new signature for fn
    new_params = sig.parameters.copy()
    new_params["to_hz"] = inspect.Parameter(
        "to_hz", inspect.Parameter.POSITIONAL_OR_KEYWORD, default=False
    )
    new_sig = inspect.Signature(
        parameters=new_params.values(), return_annotation=sig.return_annotation
    )
    fn.__signature__ = new_sig

    @preserve_signature
    @functools.wraps(fn)
    def wrapper(*args, to_hz=False, **kwargs):
        _result = fn(*args, **kwargs)
        if to_hz:
            return _result.to(u.Hz, equivalencies=[(u.cy / u.s, u.Hz)])
        return _result

    added_doc_bit = """
    Other Parameters
    ----------------
    to_hz: bool
        Set `True` to to convert function output from angular frequency to Hz
    """
    if wrapper.__doc__ is not None:
        wrapper.__doc__ += added_doc_bit
    else:
        wrapper.__doc__ = added_doc_bit

    return wrapper

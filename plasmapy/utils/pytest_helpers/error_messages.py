from astropy import units as u
from typing import Callable, Any, Dict

__all__ = [
    "call_string",
]


def _format_quantity(arg):
    formatted = f'{arg.value}'
    for base, power in zip(arg.unit.bases, arg.unit.powers):
        if power == -1:
            formatted += f"/u.{base}"
        elif power == 1:
            formatted += f"*u.{base}"
        else:
            formatted += f"*u.{base}**{power}"
    return formatted


def _format_arg(arg):
    if hasattr(arg, '__name__'):
        return arg.__name__
    elif isinstance(arg, u.quantity.Quantity):
        return _format_quantity(arg)
    else:
        return repr(arg)


def _format_kw(keyword):
    if isinstance(keyword, str):
        return str(keyword)
    elif hasattr(keyword, '__name__'):
        return keyword.__name__
    else:
        return repr(keyword)


def _exc_str(ex: Exception) -> str:
    """
    Return a string with an indefinite article and the name of
    exception ``ex``.
    """
    exception_name = ex.__name__
    use_an = exception_name[0] in 'aeiouAEIOU' and exception_name[0:4] != "User"
    indefinite_article = 'an' if use_an else 'a'
    return f"{indefinite_article} {exception_name}"


def _represent_result(result: Any) -> str:
    return f"{result.__name__}" if hasattr(result, "__name__") else f"{repr(result)}"


def call_string(f: Callable, args: Any = tuple(), kwargs: Dict = {}) -> str:
    """Return a string with the equivalent call of a function."""

    args = args if isinstance(args, tuple) else (args,)

    args_and_kwargs = ""

    for arg in args:
        args_and_kwargs += f"{_format_arg(arg)}, "

    for kwarg in kwargs:
        args_and_kwargs += f"{_format_kw(kwarg)}={_format_arg(kwargs[kwarg])}, "

    if args_and_kwargs[-2:] == ", ":
        args_and_kwargs = args_and_kwargs[:-2]

    return f"{f.__name__}({args_and_kwargs})"

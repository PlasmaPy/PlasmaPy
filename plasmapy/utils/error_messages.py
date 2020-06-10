from typing import Any, Callable, Dict

import colorama
from astropy import units as u

__all__ = ["call_string"]

_bold = colorama.Style.BRIGHT
_magenta = colorama.Fore.MAGENTA
_blue = colorama.Fore.BLUE
_cyan = colorama.Fore.CYAN
_red = colorama.Fore.RED

_exception_color = f"{_magenta}{_bold}"
_type_color = f"{_magenta}{_bold}"
_func_color = f"{_cyan}{_bold}"
_result_color = f"{_blue}{_bold}"
_message_color = f"{_red}{_bold}"


def _format_quantity(arg):
    formatted = f"{arg.value}"
    for base, power in zip(arg.unit.bases, arg.unit.powers):
        if power == -1:
            formatted += f"/u.{base}"
        elif power == 1:
            formatted += f"*u.{base}"
        else:
            formatted += f"*u.{base}**{power}"
    return formatted


def _format_arg(arg):
    if hasattr(arg, "__name__"):
        return arg.__name__
    elif isinstance(arg, u.quantity.Quantity):
        return _format_quantity(arg)
    else:
        return repr(arg)


def _format_kw(keyword):
    if isinstance(keyword, str):
        return str(keyword)
    elif hasattr(keyword, "__name__"):
        return keyword.__name__
    else:
        return repr(keyword)


def _exc_str(ex: Exception, color=_exception_color) -> str:
    """
    Return a string with an indefinite article and the name of
    exception `ex`.
    """
    if color is None:
        color = ""
        return_color = ""
    else:
        return_color = _message_color
    exception_name = ex.__name__
    use_an = exception_name[0] in "aeiouAEIOU" and exception_name[0:4] != "User"
    article = "an" if use_an else "a"
    return f"{article} {color}{exception_name}{return_color}"


def _represent_result(result: Any, color=_result_color) -> str:
    if color is None:
        color = ""
        return_color = ""
    else:
        return_color = _message_color

    if hasattr(result, "__name__"):
        return f"{color}{result.__name__}{return_color}"
    else:
        return f"{color}{repr(result)}{return_color}"


def call_string(
    f: Callable, args: Any = tuple(), kwargs: Dict = {}, color="", return_color=""
) -> str:
    """Return a string with the equivalent call of a function."""

    if color and not return_color:
        return_color = _message_color

    if not isinstance(args, tuple):
        args = (args,)

    args_and_kwargs = ""

    for arg in args:
        args_and_kwargs += f"{_format_arg(arg)}, "

    for kwarg in kwargs:
        args_and_kwargs += f"{_format_kw(kwarg)}={_format_arg(kwargs[kwarg])}, "

    if args_and_kwargs[-2:] == ", ":
        args_and_kwargs = args_and_kwargs[:-2]

    return f"{color}{f.__name__}({args_and_kwargs}){return_color}"

import pytest
from typing import Callable, Dict, Tuple, Any
import inspect
''

def _function_call_string(f: Callable, args: Tuple, kwargs: Dict) -> str:
    ...


def run_test_of_function(f: Callable, args: Any, kwargs: Dict, expected: Any):
    """Tests a function or class to make sure it returns the expected result,
    raises the appropriate exception, or issues an appropriate warning.

    Parameters
    ----------

    f: `function` or `class`
        The callable to be tested.

    args: `tuple` or `object`
        The positional arguments to `f`.

    kwargs: `dict`
        The keyword arguments to `f`.

    expected: `object`
        The expected result, exception, or warning from `f(*args, **kwargs)`.
        This may also be a tuple of length two
    """

    if not isinstance(args, tuple):
        args = (args)

    call_string = _function_call_string(f, args, kwargs)

    if inspect.isclass(expected) and issubclass(expected, Exception):

        article = 'a'

        missing_exception_errmsg = (
            f"When testing {f.__module__}, the command\n\n"
            f"    {call_string}\n\n"
            f"did not raise {article} {expected.__name__} as expected. "
        )

        with pytest.raises(expected, message=(missing_exception_errmsg)):
            f(*args, **kwargs)

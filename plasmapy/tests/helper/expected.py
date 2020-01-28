import inspect

__all__ = ["ExpectedTestOutcome"]


def _is_warning(obj) -> bool:
    """Return `True` if the argument is a warning, and `False` otherwise."""

    return inspect.isclass(obj) and issubclass(obj, Warning)


def _is_exception(obj) -> bool:
    """Return `True` if the argument is an exception, and `False` otherwise."""

    if not inspect.isclass(obj):
        return False
    else:
        return issubclass(obj, BaseException) and not issubclass(obj, Warning)


def _is_warning_and_value(obj) -> bool:
    """
    Return `True` if the argument is a `tuple` or `list` containing two
    items: a warning and an `object` that is not a warning; and `False`
    otherwise.
    """

    if not isinstance(obj, (list, tuple)) or len(obj) != 2:
        return False
    return _is_warning(obj[0]) ^ _is_warning(obj[1])


class ExpectedTestOutcome:
    """
    A class to represent the outcome of a test.

    Parameters
    ----------
    expected
        The value that is expected to be returned, the exception that
        is expected to be raised, the warning that is expected to be
        issued, or a `tuple` or `list` containing a warning and the
        expected value (in either order).

    Examples
    --------
    The most common use of this test is to represent the returned value.

    >>> outcome_is_value = ExpectedTestOutcome(42)

    Occasionally, we will want to check that a warning will be issued
    or an exception will be raised.

    >>> outcome_is_warning = ExpectedTestOutcome(RuntimeWarning)
    >>> outcome_is_exception = ExpectedTestOutcome(ValueError)

    Instances of this class contain boolean attributes to check on what
    sort of outcome is expected.

    >>> outcome_is_value.expecting_a_value
    True
    >>> outcome_is_warning.expecting_an_exception
    False
    >>> outcome_is_exception.expecting_a_warning
    False

    The expected outcome can be accessed in the ``expected_value``,
    ``expected_warning``, and ``expected_exception`` attributes.

    >>> outcome_is_value.expected_value
    42
    >>> outcome_is_warning.expected_warning
    <class 'RuntimeWarning'>
    >>> outcome_is_exception.expected_exception
    <class 'ValueError'>

    Occasionally, we will want to check that a test issues a warning
    while returning a particular value.  This can be done by passing
    a `tuple` or `list` containing a warning and the expected
    outcome (in either order).

    >>> outcome_is_value_and_warning = ExpectedTestOutcome((UserWarning, 87))
    >>> outcome_is_value_and_warning.expected_value
    87
    >>> outcome_is_value_and_warning.expected_warning
    <class 'UserWarning'>
    """

    def __init__(self, expected):

        self.expected_outcome = expected

    @property
    def expected_outcome(self):
        """
        The expected outcome of the test, which can be an exception,
        a warning, the resulting object, or a tuple that contains
        a warning and the resulting object.
        """

        if self.expecting_an_exception:
            return self.expected_exception
        elif self.expecting_a_warning and not self.expecting_a_value:
            return self.expected_warning
        elif self.expecting_a_warning and self.expecting_a_value:
            return self.expected_warning, self.expected_value
        else:
            return self.expected_value

    @expected_outcome.setter
    def expected_outcome(self, expected):

        self._info = dict()
        if _is_warning(expected):
            self._info["warning"] = expected
        elif _is_exception(expected):
            self._info["exception"] = expected
        elif _is_warning_and_value(expected):
            warning_is_first = _is_warning(expected[0])
            warning_index, value_index = (0, 1) if warning_is_first else (1, 0)
            self._info["warning"] = expected[warning_index]
            self._info["value"] = expected[value_index]
        else:
            self._info["value"] = expected

    @property
    def expecting_a_value(self) -> bool:
        """
        Return `True` if the test should return a value, and `False` otherwise.
        """

        return "value" in self._info.keys()

    @property
    def expected_value(self):
        """
        If the test is expected to return a value, then return the
        expected value.  Otherwise, raise a `RuntimeError`.
        """

        if self.expecting_a_value:
            return self._info["value"]
        else:
            raise RuntimeError("The test is not expected to return a value.")

    @property
    def expecting_an_exception(self) -> bool:
        """
        Return `True` if the test should raise an exception, and `False` otherwise.
        """

        return "exception" in self._info.keys()

    @property
    def expected_exception(self) -> BaseException:
        """
        If an exception is expected to be raised, then return that
        exception. Otherwise, raise a `RuntimeError`.
        """

        if self.expecting_an_exception:
            return self._info["exception"]
        else:
            raise RuntimeError("The test is not expected to raise an exception.")

    @property
    def expecting_a_warning(self) -> bool:
        """
        Return `True` if the test should issue a warning, and `False` otherwise.
        """

        return "warning" in self._info.keys()

    @property
    def expected_warning(self) -> Warning:
        """
        If the test is expected to issue a warning, then return that warning.
        Otherwise, raise a `RuntimeError`.
        """

        if self.expecting_a_warning:
            return self._info["warning"]
        else:
            raise RuntimeError("The test is not expected to issue a warning.")

    def __repr__(self):

        return f"ExpectedTestOutcome({self.expected_outcome})"

    def __str__(self):

        return self.__repr__()

    def __len__(self):
        """
        Return the length of the expected value, if a value is an expected
        outcome.  If ``__len__`` is undefined in the expected value, then
        return ``1``.
        """

        if self.expecting_a_value:
            has_a_len = hasattr(self.expected_value, "__len__")
            return len(self.expected_value) if has_a_len else 1
        else:
            return NotImplemented

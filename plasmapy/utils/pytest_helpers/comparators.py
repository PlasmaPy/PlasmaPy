"""Contains test comparators: tools that automate comparisons."""

from typing import List
from plasmapy.utils.pytest_helpers.expected import ExpectedTestOutcome
from plasmapy.utils.pytest_helpers.actual import ActualTestOutcome


class CompareActualExpected:
    def __init__(self, actual: ActualTestOutcome, expected: ExpectedTestOutcome):
        if not isinstance(actual, ActualTestOutcome):
            raise TypeError("Expecting an instance of ActualTestOutcome")
        if not isinstance(expected, ExpectedTestOutcome):
            raise TypeError("Expecting an instance of ExpectedTestOutcome")

        self._actual = actual
        self._expected = expected

        self._make_comparison()

    @property
    def actual(self) -> ActualTestOutcome:
        """
        The instance of `~plasmapy.utils.pytest_helpers.actual.ActualTestOutcome`
        that is being compared.
        """
        return self._actual

    @property
    def expected(self) -> ExpectedTestOutcome:
        """
        The instance of `~plasmapy.utils.pytest_helpers.expected.ExpectedTestOutcome`
        that is being compared.
        """
        return self._expected

    @property
    def test_passed(self) -> bool:
        """
        Returns `True` if the actual outcome matches the expected outcome,
        and `False` otherwise.
        """
        return bool(self._error_message)

    @property
    def error_messages(self) -> List[str]:
        """
        Returns a list of the error messages
        """
        return self._error_message

    def _make_comparison(self):
        ...

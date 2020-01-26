import pytest
import collections

from plasmapy.utils.pytest_helpers.runners import (
    function_test_runner,
    attr_test_runner,
    method_test_runner,
)

from plasmapy.utils.pytest_helpers.tests.sample_functions import (
    return_42,
    return_42_meters,
    return_np_array,
    return_sum_of_two_args_and_kwargs,
    return_none,
    raise_sample_exception,
    issue_sample_warning_and_return_42,
    SampleClass,
    SampleException,
    SampleExceptionSubclass,
    SampleWarning,
    SampleWarningSubclass,
)


class FunctionTestCase:
    def __init__(
        self, expected, function, args=None, kwargs=None, exception=None, errmsg=None
    ):
        """
        Store information for a test of ``function_test_runner``.

        If no exception is provided, then the test is assumed to pass.

        """
        self.expected = expected
        self.function = function
        self.args = args if args is not None else ()
        self.kwargs = kwargs if kwargs is not None else {}
        self.expected_exception = exception
        self.test_should_pass = exception is None
        self.errmsg = "" if errmsg is None else errmsg


function_test_cases = [
    FunctionTestCase(expected=42, function=return_42),
    FunctionTestCase(expected=43, function=return_42, exception=pytest.fail.Exception),
]


@pytest.mark.parametrize("case", function_test_cases)
def test_function_test_runner(case):

    if case.test_should_pass:
        function_test_runner(case.expected, case.function, case.args, case.kwargs)
    else:
        with pytest.raises(case.expected_exception, match=case.errmsg):
            function_test_runner(case.expected, case.function, case.args, case.kwargs)

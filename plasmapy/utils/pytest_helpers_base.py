from astropy import units as u
from numbers import Real
import inspect


class Expected:
    """



    """


    def __init__(self, expected_outcome):

        self.expected_outcome = expected_outcome


    @property
    def excepted_outcome(self):
        return self._expected_outcome

    @expected_outcome.setter
    def expected_outcome(self, value):




    @property
    def exception(self):
        ...

    @exception.setter
    def exception(self, exception):



class Result:
    def __init__(self):
        ...

    def __eq__(self, other):
        ...





class RunTestBase:

    def



    @property
    def result(self):
        ...

    @property
    def expected(self, *values):
        ...


    @property
    def atol(self):
        return self._atol

    @atol.setter
    def atol(self, value: Real):


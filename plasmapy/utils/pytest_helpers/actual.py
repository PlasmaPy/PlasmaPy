import pytest
import warnings
import copy

__all__ = ["ActualTestOutcome"]


class ActualTestOutcome:
    def __init__(self, inputs):
        self.inputs = inputs
        try:
            with pytest.warns(Warning) as warnings_record:
                self.result = self.inputs.call()
                self.warnings_record = copy.copy(warnings_record)
                self.warning_was_issued = bool(self.warnings_record)
                warnings.warn(Warning)
        except Exception as exception_info:
            self.exception_was_raised = True
            self.warning_was_issued = False
            self.exception_info = copy.copy(exception_info)
        else:
            self.exception_was_raised = False
            self.exception_info = None

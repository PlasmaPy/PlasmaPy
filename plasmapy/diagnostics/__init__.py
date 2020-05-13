"""The diagnostics subpackage contains tools for experimental research.
Currently, we have functionality for analyzing data from Langmuir probes.
"""
__all__ = ['AbstractProbe', 'XAbstractDiagnostic']

import abc
import xarray as xr


class AbstractProbe(abc.ABC):
    """Abstract class for defining probe characteristics/parameters."""


class XAbstractDiagnostic(abc.ABC):
    """
    Abstract class for `xarray` flavored diagnostic classes.
    """
    def __init__(self, dataset: xr.Dataset):
        self._ds = dataset
        self._probe_class = AbstractProbe

    @property
    def is_probe_defined(self) -> bool:
        return 'probe_parameters' in self._ds.attrs.keys()

    @property
    def is_probe_correct(self) -> bool:
        return isinstance(self.probe_parameters, self._probe_class)

    @property
    def probe_parameters(self):
        try:
            return self._ds.attrs['probe_parameters']
        except KeyError:
            raise AttributeError(f"Dataset attribute 'probe_parameters' is not "
                                 f"defined.")

    @probe_parameters.setter
    def probe_parameters(self, value):
        if isinstance(value, self._probe_class):
            self._ds.attrs['probe_parameters'] = value
        else:
            raise ValueError(f"Expected instance of {self._probe_class} "
                             f"and got {value.__class__}.")

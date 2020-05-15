__all__ = [
    "AbstractProbe",
    "XAbstractDiagnostic",
    "XDiagnostics",
]

import abc
import importlib
import xarray as xr

from warnings import warn


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


# @xr.register_dataset_accessor('plasmapy')
# class PlasmaPyAccessor:
class XDiagnostics:
    __available_diagnostics = {
        'swept_langmuir': (".swept_langmuir.xdiagnostic",
                           "XSweptLangmuirDiagnostic"),
    }

    # def __init__(self, xr_obj: xr.Dataset):
    #     self._ds = xr_obj

    def __repr__(self) -> str:
        summary = [
            super().__repr__(),
            "",
            "Enabled   Available Diagnostic",
        ]

        for diag in self.__available_diagnostics.keys():
            sum_str = "    ["
            sum_str += "x" if hasattr(xr.Dataset, diag) else " "
            sum_str += f"]   {diag}"
            summary.append(sum_str)

        return '\n'.join(summary)

    @property
    def available(self):
        return list(self.__available_diagnostics)

    def enable(self, *args):
        for arg in args:
            if arg in self.__available_diagnostics:
                if hasattr(xr.Dataset, arg):
                    # Do NOT register accessor to xarray.Dataset if already there
                    continue

                diag_name = arg
                mod_name = self.__available_diagnostics[diag_name][0]
                cls_name = self.__available_diagnostics[diag_name][1]
                cls = getattr(
                    importlib.import_module(mod_name, package="plasmapy.diagnostics"),
                    cls_name,
                )
                xr.register_dataset_accessor(diag_name)(cls)
            else:
                warn(f"Requested diagnostic '{arg}' is not among available "
                     f"diagnostics.")

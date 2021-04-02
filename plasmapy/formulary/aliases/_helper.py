import importlib
import inspect
import sys


def _generate_aliases(_name_):
    module = importlib.import_module(_name_.replace(".aliases", ""))
    _all_ = []
    for name, function in inspect.getmembers(module):
        if not name.startswith("_") and name.endswith("_"):
            _all_.append(name)
            setattr(sys.modules[_name_], name, function)
    _all_.sort()
    return _all_

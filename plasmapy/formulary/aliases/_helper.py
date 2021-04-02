import functools
import importlib
import inspect
import sys


def _generate_namespace(_name_, func_condition, mod_name_replacer):
    module = importlib.import_module(mod_name_replacer(_name_))
    _all_ = []
    for name, function in inspect.getmembers(module):
        if func_condition(name):
            _all_.append(name)
            setattr(sys.modules[_name_], name, function)
    _all_.sort()
    return _all_


_generate_aliases = functools.partial(
    _generate_namespace,
    func_condition=lambda name: not name.startswith("_") and name.endswith("_"),
    mod_name_replacer=lambda modname: modname.replace(".aliases", ""),
)

_generate_lite_versions = functools.partial(
    _generate_namespace,
    func_condition=lambda name: name.endswith("lite"),
    mod_name_replacer=lambda modname: modname.replace(".lite", ""),
)

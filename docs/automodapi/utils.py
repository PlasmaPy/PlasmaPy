__all__ = ["automod_groupings", "find_mod_objs"]

import inspect
import os

from importlib import import_module
from sphinx.application import Sphinx

package_dir = os.path.abspath(os.path.dirname(__file__))
templates_dir = os.path.join(package_dir, "templates")

def automod_groupings(app: Sphinx):
    default_groupings = {
        "modules",
        "classes",
        "exceptions",
        "warnings",
        "functions",
        "variables",
    }
    custom_groups = set(app.config.automod_custom_groups)

def get_custom_grouping_info(app: Sphinx):
    try:
        _info = app.config.automod_custom_groups
    except AttributeError:
        _info = {}

    return _info


def find_mod_objs(modname: str, app: Sphinx = None):
    if app is not None:
        if isinstance(app, Sphinx):
            cgroups_def = app.config.automod_custom_groups
        else:
            # assuming dict for testing
            cgroups_def = app

        cgroups = set(cgroups_def)
    else:
        cgroups_def = {}
        cgroups = set()

    mod = import_module(modname)
    pkg_name = modname.split(".")[0]

    # define what to search
    pkg_names = {name for name in mod.__dict__.keys() if not name.startswith("_")}
    if hasattr(mod, "__all__"):
        no_all = False
        names_to_search = set(mod.__all__)
    else:
        no_all = True
        names_to_search = pkg_names

    # filter pkg_names
    for name in pkg_names.copy():
        obj = getattr(mod, name)

        if not no_all and name in names_to_search:
            continue

        ismod = inspect.ismodule(obj)
        ispkg = ismod and obj.__package__ == obj.__name__

        # remove test folders
        if ispkg and obj.__package__.split(".")[-1] == "tests":
            pkg_names.remove(name)
            continue

        # remove 3rd party objects
        if ismod and obj.__package__.split(".")[0] != pkg_name:
            pkg_names.remove(name)
            continue
        elif (
            not ismod
            and hasattr(obj, "__module__")
            and obj.__module__.split(".")[0] != pkg_name
        ):
            # Note: this will miss ufuncs like numpy.sqrt since they do not have
            #       a __module__ property
            pkg_names.remove(name)
            continue

        # remove non direct sub-pkgs and mods of modname
        if ismod:
            if not obj.__name__.startswith(modname):
                pkg_names.remove(name)
                continue
            else:
                nm = obj.__name__[len(modname):].split(".")
                nm.remove("")
                if len(nm) != 1:
                    pkg_names.remove(name)
                    continue

    # find local modules first
    names_of_modules = set()
    for name in pkg_names.copy():
        obj = getattr(mod, name)

        if inspect.ismodule(obj):
            names_of_modules.add(name)

    mod_objs = {"modules": {"names": []}}
    if len(names_of_modules) > 0:
        names_of_modules = names_of_modules
        mod_objs["modules"]["names"] = list(names_of_modules)
        names_to_search = names_to_search - names_of_modules

    # find and filter custom groups
    for name in cgroups:
        dunder = cgroups_def[name]["dunder"]
        if hasattr(mod, dunder):
            custom_names = set(getattr(mod, dunder))
        else:
            continue

        if len(custom_names) > 0:
            mod_objs.update({name: {"names": list(custom_names)}})
            names_to_search = names_to_search - custom_names

    # gather all remaining groups
    mod_objs.update(
        {
            "classes": {"names": []},
            "exceptions": {"names": []},
            "warnings": {"names": []},
            "functions": {"names": []},
            "variables": {"names": []},
        }
    )
    for name in names_to_search:
        obj = getattr(mod, name)

        if inspect.isroutine(obj):
            # is a user-defined or built-in function
            mod_objs["functions"]["names"].append(name)
        elif inspect.isclass(obj):
            if issubclass(obj, Warning):
                mod_objs["warnings"]["names"].append(name)
            elif issubclass(obj, BaseException):
                mod_objs["exceptions"]["names"].append(name)
            else:
                mod_objs["classes"]["names"].append(name)
        else:
            mod_objs["variables"]["names"].append(name)

    # retrieve and defined qualnames and objs
    for obj_type in list(mod_objs):
        if len(mod_objs[obj_type]["names"]) == 0:
            del mod_objs[obj_type]
            continue

        mod_objs[obj_type].update(
            {
               "qualnames": [],
               "objs": [],
            },
        )
        for name in list(mod_objs[obj_type]["names"]):
            # Note:  The 'qualname' is always constructed with 'name' so when
            #        something like
            #
            #        def func(...):
            #            ...
            #
            #        f2 = func
            #
            #        is done, then the 'qualname' ends with 'f2' and not 'func'.
            #
            obj = getattr(mod, name)

            ismod = inspect.ismodule(obj)
            ispkg = ismod and obj.__package__ == obj.__name__

            if not ismod and no_all:
                # only allow local objects to be collected
                # - at this point modules & pkgs should have already been
                #   filtered for direct sub-modules and pkgs
                if not hasattr(obj, "__module__"):
                    # this would be a locally defined variable like
                    # plasmapy.__citation__
                    pass
                elif obj.__module__ != modname:
                    mod_objs[obj_type]["names"].remove(name)
                    continue

            if ismod:
                obj_renamed = obj.__name__.split(".")[-1] != name
            elif not hasattr(obj, "__name__"):
                obj_renamed = False
            else:
                obj_renamed = obj.__name__ != name

            if ismod and obj_renamed:
                qualname = f"{obj.__package__}.{name}"
            elif ismod and not obj_renamed:
                qualname = obj.__name__
            elif obj_renamed or not hasattr(obj, "__module__"):
                # can not tell if the object was renamed in modname or in
                # obj.__module__, so assumed it happened in modname
                qualname = f"{modname}.{name}"
            elif obj.__module__.split(".")[0] != pkg_name:
                # this will catch scenarios like typing alias definitions where
                # __module__ == typing even when defined locally
                qualname = f"{modname}.{name}"
            else:
                qualname = f"{obj.__module__}.{name}"

            mod_objs[obj_type]["qualnames"].append(qualname)
            mod_objs[obj_type]["objs"].append(obj)

    return mod_objs

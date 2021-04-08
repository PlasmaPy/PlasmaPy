__all__ = ["find_mod_objs"]

import inspect
import sys

from sphinx.application import Sphinx


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

    return default_groupings, custom_groups


def find_mod_objs(modname: str, only_locals=False, app: Sphinx = None):
    if app is not None:
        if isinstance(app, Sphinx):
            cgroups_def = app.config.automod_custom_groups
            # dgroups, cgroups = automodapi_groupings(app)
        else:
            # assuming dict for testing
            cgroups_def = app

        cgroups = set(cgroups_def)
    else:
        cgroups_def = {}
        cgroups = set()

    mod_objs = {}

    __import__(modname)
    mod = sys.modules[modname]

    # define what to search
    pkg_names = {name for name in mod.__dict__.keys() if not name.startswith("_")}
    if hasattr(mod, "__all__"):
        all_names = set(mod.__all__)
    else:
        all_names = pkg_names
    names_to_search = all_names

    # find local modules first
    names_of_modules = []
    for name in pkg_names:
        obj = getattr(mod, name)

        if inspect.ismodule(obj) and obj.__spec__.parent == modname:
            names_of_modules.append(name)
    mod_objs.update({"modules": {"names": []}})
    if len(names_of_modules) > 0:
        names_of_modules = set(names_of_modules)
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
            if inspect.ismodule(obj):
                if obj.__package__ == obj.__name__:
                    # this is a directory
                    qualname = ".".join(obj.__package__.split(".")[:-1] + [name])
                else:
                    # this is a py file
                    qualname = f"{obj.__package__}.{name}"

            elif hasattr(obj, "__module__"):
                qualname = f"{obj.__module__}.{name}"
            else:
                qualname = f"{modname}.{name}"

            if only_locals:
                allowed = modname
                if isinstance(only_locals, (tuple, list)):
                    allowed = tuple(only_locals)
                if not qualname.startswith(allowed):
                    del mod_objs[obj_type]["names"][name]
                    continue

            mod_objs[obj_type]["qualnames"].append(qualname)
            mod_objs[obj_type]["objs"].append(obj)

    return mod_objs

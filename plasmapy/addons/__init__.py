__all__ = []

from pkg_resources import (iter_entry_points, resource_filename)


def _parse_readme():
    fname = resource_filename('plasmapy', 'addons/README.md')
    with open(fname, encoding='utf-8') as f:
        readme = f.read()

    return readme


__readme__ = _parse_readme()


# import addon entry points
for ep in iter_entry_points('plasmapy.addons'):
    # __dict__[ep.name] = ep.load()
    globals()[ep.name] = ep.load()
    __all__.append(ep.name)

del ep, iter_entry_points, resource_filename

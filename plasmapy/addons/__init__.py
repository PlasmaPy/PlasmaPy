

from pkg_resources import resource_filename


def _parse_readme():
    fname = resource_filename('plasmapy', 'addons/README.md')
    with open(fname, encoding='utf-8') as f:
        readme = f.read()

    return readme


def about():
    print(_parse_readme())

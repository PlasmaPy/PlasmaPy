from pkg_resources import resource_filename


def about():
    fname = resource_filename('plasmapy', 'addons/README.md')
    with open(fname, encoding='utf-8') as f:
        print(f.read())

import astropy.units as u

class Example:
    def __init__(self, x):
        self._x = x

    @property
    def x(self):
        return u.Quantity(self._x, u.m)

    @x.setter
    def x(self, x):
        if isinstance(x, u.Quantity):
            self._x = x.si.value
        else:
            self._x = x

e = Example(0)
assert e.x == 0 * u.m
assert e._x == 0
e._x = 2
assert e.x == 2 * u.m
assert e._x == 2
e.x = 3
assert e.x == (3 * u.m)
assert e._x == 3

def AstropyBackedField(self, variable):

    return None

class AstropyBackedExample:
    def __init__(self, x: u.m):
        self.x = AstropyBackedField(self, x)

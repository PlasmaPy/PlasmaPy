"""Defines the UncertaintyQuantity class."""

import astropy.units as u
import numpy as np

class UncertaintyQuantity(u.Quantity):

    def set_uncertainty(self, uncertainty):
        self.uncertainty = uncertainty
        return self

    def relative(self):
        return float(self.uncertainty / u.Quantity(self))

    def __new__(self, base, uncertainty):
        new = u.Quantity(base)
        new.__class__ = UncertaintyQuantity
        new.uncertainty = uncertainty
        return new

    def base(self):
        return self.get_base()

    def get_base(self):
        result = u.Quantity(self)
        return result

    def get_uncertainty(self):
        return u.Quantity(self.uncertainty)

    def to(self, *vars, **kwargs):
        result = super(UncertaintyQuantity, self).to(*vars, **kwargs)
        result.uncertainty = self.uncertainty.to(*vars, **kwargs)
        return result

    def __eq__(self, other):
        if not u.Quantity(self) == u.Quantity(other):
            return False

        if not u.Quantity(self.uncertainty) == u.Quantity(other.uncertainty):
            return False

        return True

    @staticmethod
    def _uncertainty_method(a, b, term1, term2, **kwargs):

        if(type(a) != UncertaintyQuantity or a.uncertainty == 0):
            return term2()

        if(type(b) != UncertaintyQuantity or b.uncertainty == 0):
            return term1()

        return np.sqrt(np.power(term1(), 2, **kwargs) + np.power(term2(), 2, **kwargs), **kwargs)

    @staticmethod
    def get_pow_uncertainty(a, b, **kwargs):

        def term1():
            return (np.abs(u.Quantity(b) *
                           np.power(u.Quantity(a), u.Quantity(b) - 1, **kwargs), **kwargs) *
                     u.Quantity(a.uncertainty))

        def term2():
            return (np.abs(np.power(u.Quantity(a), u.Quantity(b), **kwargs) *
                           np.log(u.Quantity(a), **kwargs), **kwargs) *
                    u.Quantity(b.uncertainty))

        return UncertaintyQuantity._uncertainty_method(a, b, term1, term2, **kwargs)

    @staticmethod
    def get_add_uncertainty(a, b, **kwargs):

        def term1():
            return u.Quantity(a.uncertainty)

        def term2():
            return u.Quantity(b.uncertainty)

        return UncertaintyQuantity._uncertainty_method(a, b, term1, term2, **kwargs)

    @staticmethod
    def get_sub_uncertainty(a, b, **kwargs):
        return UncertaintyQuantity.get_add_uncertainty(a, b, **kwargs)

    @staticmethod
    def get_mul_uncertainty(a, b, **kwargs):

        def term1():
            return np.abs(u.Quantity(b)) * u.Quantity(a.uncertainty)

        def term2():
            return np.abs(u.Quantity(a)) * u.Quantity(b.uncertainty)

        return UncertaintyQuantity._uncertainty_method(a, b, term1, term2, **kwargs)

    @staticmethod
    def get_truediv_uncertainty(a, b, **kwargs):

        def term1():
            return (np.abs(1 / u.Quantity(b), **kwargs) *
                    u.Quantity(a.uncertainty))

        def term2():
            return (np.abs(u.Quantity(a) / u.Quantity(b) ** 2, **kwargs) *
                     u.Quantity(b.uncertainty))

        return UncertaintyQuantity._uncertainty_method(a, b, term1, term2, **kwargs)

    def __truediv__(self, other):
        result = super(UncertaintyQuantity, self).__truediv__(other)

        if not (type(other) == UncertaintyQuantity or type(other) == u.Quantity):
            other = u.Quantity(1 * other)

        result.uncertainty = UncertaintyQuantity.get_truediv_uncertainty(self, other)

        return result

    def __pow__(self, other):
        result = super(UncertaintyQuantity, self).__pow__(other)

        if not (type(other) == UncertaintyQuantity or type(other) == u.Quantity):
            other = u.Quantity(1 * other)

        result.uncertainty = UncertaintyQuantity.get_pow_uncertainty(self, other)

        return result

    def __add__(self, other):
        result = super(UncertaintyQuantity, self).__add__(other)

        if not (type(other) == UncertaintyQuantity or type(other) == u.Quantity):
            other = u.Quantity(1 * other)

        result.uncertainty = UncertaintyQuantity.get_add_uncertainty(self, other)

        return result

    def __sub__(self, other):
        result = super(UncertaintyQuantity, self).__sub__(other)

        if not (type(other) == UncertaintyQuantity or type(other) == u.Quantity):
            other = u.Quantity(1 * other)

        result.uncertainty = UncertaintyQuantity.get_sub_uncertainty(self, other)

        return result

    def __mul__(self, other):
        result = super(UncertaintyQuantity, self).__mul__(other)

        if not (type(other) == UncertaintyQuantity or type(other) == u.Quantity):
            other = u.Quantity(1 * other)

        result.uncertainty = UncertaintyQuantity.get_mul_uncertainty(self, other)

        return result

    def __array_ufunc__(self, *vars, **kwargs):
        result = super(UncertaintyQuantity, self).__array_ufunc__(*vars, **kwargs)

        if(vars[0].__name__ == 'sqrt'):
            this = vars[2]
            result.uncertainty = UncertaintyQuantity.get_pow_uncertainty(this, 0.5, **kwargs)

        if(vars[0].__name__ == 'power'):
            this = vars[2]
            other = vars[3]
            result.uncertainty = UncertaintyQuantity.get_pow_uncertainty(this, other, **kwargs)

        if(vars[0].__name__ == 'add'):
            this = vars[2]
            other = vars[3]
            result.uncertainty = UncertaintyQuantity.get_add_uncertainty(this, other, **kwargs)

        if(vars[0].__name__ == 'subtract'):
            this = vars[2]
            other = vars[3]
            result.uncertainty = UncertaintyQuantity.get_sub_uncertainty(this, other, **kwargs)

        if(vars[0].__name__ == 'multiply'):
            this = vars[2]
            other = vars[3]
            result.uncertainty = UncertaintyQuantity.get_mul_uncertainty(this, other, **kwargs)

        if(vars[0].__name__ == 'true_divide'):
            this = vars[2]
            other = vars[3]
            result.uncertainty = UncertaintyQuantity.get_truediv_uncertainty(this, other, **kwargs)

        if(vars[0].__name__ == 'square'):
            this = vars[2]
            result.uncertainty = UncertaintyQuantity.get_pow_uncertainty(this, 2, **kwargs)

        return result

    def __str__(self):
        return super(UncertaintyQuantity, self).value.__str__() + \
            " ± " + \
            self.uncertainty.__str__()
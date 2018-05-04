"""Defines the UncertaintyQuantity class."""

import astropy.units as u
import numpy as np


class UncertaintyQuantity(u.Quantity):
    """
    Extension of the Quantity class to implement uncertainty properties and propagation.

    Attributes
    ----------
    uncertainty : `astropy.units.Quantity`
        The uncertainty in the value of the parent object. This must be a Quantity to prevent
        recursion.

    Notes
    -----
    UncertaintyQuantity is highly experimental. Many but the most common operators and numpy
    functions are not supported. Currently supported operations are:

    - addition
    - subtraction
    - multiplication
    - division
    - power
    - sqrt
    - square

    """

    def __new__(self, base, uncertainty):
        new = u.Quantity(base)
        new.__class__ = UncertaintyQuantity
        new._uncertainty = u.Quantity(uncertainty)
        return new

    @property
    def uncertainty(self):
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):

        if(value < 0):
            raise ValueError('Uncertainty cannot be negative.')

        self._uncertainty = value.to(self.unit)

    def relative(self):
        return float(self.uncertainty / u.Quantity(self))

    def base(self):
        return self.get_base()

    def get_base(self):
        result = u.Quantity(self)
        return result

    def to(self, *vars, **kwargs):
        result = super(UncertaintyQuantity, self).to(*vars, **kwargs)
        result._uncertainty = self.uncertainty.to(*vars, **kwargs)
        return result

    def __eq__(self, other):
        """
        Equality is defined based on both the parent Quantity and the uncertainty Quantity.

        """

        if not u.Quantity(self) == u.Quantity(other):
            return False

        if not u.Quantity(self._uncertainty) == u.Quantity(other._uncertainty):
            return False

        return True

    def __truediv__(self, other):
        """Redirects to the numpy operator"""
        return np.true_divide(self, other)

    def __pow__(self, other):
        """Redirects to the numpy operator"""
        return np.power(self, other)

    def __add__(self, other):
        """Redirects to the numpy operator"""
        return np.add(self, other)

    def __sub__(self, other):
        """Redirects to the numpy operator"""
        return np.subtract(self, other)

    def __mul__(self, other):
        """Redirects to the numpy operator"""
        return np.multiply(self, other)

    derivative_dict = {
            np.add:         [lambda a, b: 1,
                             lambda a, b: 1],
            np.subtract:    [lambda a, b: 1,
                             lambda a, b: 1],
            np.multiply:    [lambda a, b: b,
                             lambda a, b: a],
            np.true_divide: [lambda a, b: 1 / b,
                             lambda a, b: a / b ** 2],
            np.square:      [lambda a: 2 * a],
            np.sqrt:        [lambda a: 0.5 / a ** 0.5],
            np.power:       [lambda a, b: b * a ** (b - 1),
                             lambda a, b: np.log(a) * a ** b]}

#    @staticmethod
#    def absolute_uncertainty_rule(terms):
#        """
#        The uncertainty propagation rule for 100% confidence intervals:
#        > Df(x, y) = abs(df/dx) * Dx + abs(df/dy) * Dy
#
#        """
#
#        return np.sum(np.abs(term) for term in terms)

    @staticmethod
    def general_uncertainty_propagation_rule(terms):
        """
        The standard uncertainty propagation rule for 68% confidence intervals:
        > Df(x, y) = sqrt(abs(df/dx)^2 * Dx^2 + abs(df/dy)^2 * Dy^2)

        """

        return np.sqrt(np.sum((term**2 for term in terms)))

    def __array_ufunc__(self, *vars, **kwargs):
        """
        Implements uncertainty propagation by extending __array_ufunc__, which is called whenever
        an ufunc is applied. If not applicable or not implemented the normal result of the ufunc
        is returned. Otherwise the standard rules of uncertainty calculus are applied to obtain the
        uncertainty of the function result.

        """

        # This line unpacks one or both parameters of the operation into the variable params
        ufunc, _, *params = vars

        # Defines a function to convert any input value to a Quantity
        def toQuantity(value):
            if not (type(value) == UncertaintyQuantity or type(value) == u.Quantity):
                return u.Quantity(1 * value)
            else:
                return u.Quantity(value)

        # Convert the input parameters to Quantity objects
        quantity_params = [toQuantity(param) for param in params]

        # Execute the ufunc to obtain the normal result of the operation
        result = ufunc(*quantity_params, **kwargs)

        # Return the result if no derivative uncertainty function is defined
        if(ufunc not in self.derivative_dict):
            return result

        # Ensure the result of the normal ufunc is an UncertaintyQuantity
        result.__class__ = UncertaintyQuantity

        # Obtain the derivative lambdas corresponding to the ufunc
        derivative_funcs = self.derivative_dict[ufunc]

        # Generate a list of uncertainty terms. The terms are defined as the product of the
        # absolute value of the derivative to a parameter and the uncertainty of that parameter.
        terms = []
        for i in np.arange(len(params)):

            # Do not calculate the term if the parameter has no uncertainty or it is zero
            if(params[i].__class__ == UncertaintyQuantity and params[i]._uncertainty != 0):

                # Add the term to the list
                terms.append(np.abs(derivative_funcs[i](*quantity_params)) * params[i].uncertainty)

        # Several methods are available to obtain the uncertainty based on the nature of the
        # interval. Needs some easy method of switching, ie. subclassing UncertaintyQuantity.
        result._uncertainty = UncertaintyQuantity.general_uncertainty_propagation_rule(terms)

        return result

    def __str__(self):
        """Generates a shortened representation of the UncertaintyQuantity object"""

        return '{0} ± {1}{2:s}'.format(
                self.value,
                self._uncertainty.value,
                self._unitstr)

    def __repr__(self):
        """Generates a string representation of the UncertaintyQuantity object"""

        prefixstr = '<' + self.__class__.__name__ + ' '

        return '{0}{1} ± {2}{3:s}>'.format(
                prefixstr,
                self.value,
                self._uncertainty.value,
                self._unitstr)

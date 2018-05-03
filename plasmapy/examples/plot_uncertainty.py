# coding: utf-8

from plasmapy.uncertainty import UncertaintyQuantity
import astropy.units as u
import numpy as np

value0 = UncertaintyQuantity(30 * u.m, 5 * u.m)
value1 = UncertaintyQuantity(10 * u.m, 1 * u.m)

print("Value a = " + str(value1))

print("Base of a = " + str(value1.base()))

print("Uncertainty of a = " + str(value1.uncertainty))

print("Value b = " + str(value1))

summation = value0 + value1
print("Summation (a + b) = " + str(summation))

subtraction = value0 - value1
print("Subtraction (a - b) = " + str(subtraction))

multiplication = value0 * value1
print("Multiplication (a * b) = " + str(multiplication))

division = value0 / value1
print("Division (a / b) = " + str(division))

power = (value0 / u.m) ** UncertaintyQuantity(2, 0.1)
print("Power (a ** (2 ± 0.1)) = " + str(power))

root = np.sqrt(value0)
print("Root (sqrt(a)) = " + str(root))

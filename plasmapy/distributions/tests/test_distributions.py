import astropy.units as u
import numpy as np

from plasmapy.distributions.distributions import *


def test_distributions_creation():
    # Smoke test distribution funciton creation
    data = np.random.rand(3, 3, 3)
    indexes = [np.random.rand(3) for i in range(3)]
    index_types = ['time', 'space', 'velocity']
    DiscreteDistributionFunction(data, indexes, index_types)

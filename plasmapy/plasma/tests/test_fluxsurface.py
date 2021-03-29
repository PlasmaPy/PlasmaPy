import numpy as np
import pytest


def test_fs_Bmag(num_regression, flux_surface):
    num_regression.check(
        {"fs-averaged mod B": flux_surface.flux_surface_average(flux_surface.Bmag)}
    )

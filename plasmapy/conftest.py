'''Adds slow marker and forces matplotlib to use non-gui backends during tests'''
import pytest

# Force MPL to use non-gui backends for testing.
try:
    import matplotlib
except ImportError:
    pass
else:
    matplotlib.use('Agg')

# Add slow marker
# coverage : ignore
def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        ("slow: mark test as slow to run. Used to mark tests that execute in more than 1000x the "
        "time of the median test.")
    )

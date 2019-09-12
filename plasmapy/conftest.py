    # Force MPL to use non-gui backends for testing.
import pytest

try:
    import matplotlib
except ImportError:
    pass
else:
    matplotlib.use('Agg')

def pytest_addoption(parser):
    parser.addoption(
        "--not-slow", action="store_true", default=False, help='Pytest will skip slow tests.'
    )
    parser.addoption(
        "--slow", action="store_true", default=False, help='Pytest will only run slow tests.'
    )

def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        ("slow: mark test as slow to run. Test can be skipped with --not-slow or exclusively "
        "executed with --slow.")
    )

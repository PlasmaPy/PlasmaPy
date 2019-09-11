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

def pytest_collection_modifyitems(config, items):
    skip_condtion = None
    if config.getoption("--not-slow") and config.getoption("--slow"):
        # User wants to run both the not-slow tests and the slow tests, which is the same as running
        #   with no options
        pass
    elif config.getoption("--not-slow"):
        # Skip slow tests
        skip_mark = pytest.mark.skip(reason="Test is marked slow.")
        skip_condition = lambda x: "slow" in x.keywords
    elif config.getoption("--slow"):
        skip_mark = pytest.mark.skip(reason="Test isn't marked slow.")
        skip_condition = lambda x: "slow" not in x.keywords
    if skip_condition is not None:
        for item in items:
            if skip_condition(item):
                item.add_marker(skip_mark)

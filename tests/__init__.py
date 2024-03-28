"""Tests for PlasmaPy."""

# The empty __init__.py files in tests/ and subdirectories make each
# directory be treated as containing packages. Removing the __init__.py
# files may lead to problems with per-file ignores specified in mypy.ini
# and problems with importing tests if there are files with the same
# name in multiple directories (like conftest.py). The __init__.py files
# are also necessary if a test file imports anything from a different
# file in the directory.

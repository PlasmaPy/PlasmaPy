"""Test functionality related to setup.py and importing plasmapy."""

from .test_import_helpers import (
    test_check_versions_newer,
    test_check_versions_older,
)

from .test_availability import (
    test_environment_txt_and_yml_files,
    test_for_no_duplicate_lines,
    test_no_formatted_string_literals,
    test_requirements_issubset_automated_code_tests,
    test_requirements_issubset_of_other_requirements_files,
)

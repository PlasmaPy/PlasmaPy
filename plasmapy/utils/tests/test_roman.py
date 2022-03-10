import numpy as np
import pytest

import plasmapy.utils.roman as roman

from plasmapy.utils.pytest_helpers import run_test

ints_and_roman_numerals = [
    (1, "I"),
    (2, "II"),
    (3, "III"),
    (4, "IV"),
    (5, "V"),
    (6, "VI"),
    (7, "VII"),
    (8, "VIII"),
    (9, "IX"),
    (10, "X"),
    (11, "XI"),
    (12, "XII"),
    (13, "XIII"),
    (14, "XIV"),
    (15, "XV"),
    (16, "XVI"),
    (17, "XVII"),
    (18, "XVIII"),
    (19, "XIX"),
    (20, "XX"),
    (21, "XXI"),
    (22, "XXII"),
    (23, "XXIII"),
    (24, "XXIV"),
    (25, "XXV"),
    (26, "XXVI"),
    (27, "XXVII"),
    (28, "XXVIII"),
    (29, "XXIX"),
    (30, "XXX"),
    (31, "XXXI"),
    (32, "XXXII"),
    (33, "XXXIII"),
    (34, "XXXIV"),
    (35, "XXXV"),
    (36, "XXXVI"),
    (37, "XXXVII"),
    (38, "XXXVIII"),
    (39, "XXXIX"),
    (40, "XL"),
    (41, "XLI"),
    (42, "XLII"),
    (43, "XLIII"),
    (44, "XLIV"),
    (45, "XLV"),
    (46, "XLVI"),
    (47, "XLVII"),
    (48, "XLVIII"),
    (49, "XLIX"),
    (50, "L"),
    (51, "LI"),
    (52, "LII"),
    (53, "LIII"),
    (54, "LIV"),
    (55, "LV"),
    (56, "LVI"),
    (57, "LVII"),
    (58, "LVIII"),
    (59, "LIX"),
    (60, "LX"),
    (61, "LXI"),
    (62, "LXII"),
    (63, "LXIII"),
    (64, "LXIV"),
    (65, "LXV"),
    (66, "LXVI"),
    (67, "LXVII"),
    (68, "LXVIII"),
    (69, "LXIX"),
    (70, "LXX"),
    (71, "LXXI"),
    (72, "LXXII"),
    (73, "LXXIII"),
    (74, "LXXIV"),
    (75, "LXXV"),
    (76, "LXXVI"),
    (77, "LXXVII"),
    (78, "LXXVIII"),
    (79, "LXXIX"),
    (80, "LXXX"),
    (81, "LXXXI"),
    (82, "LXXXII"),
    (83, "LXXXIII"),
    (84, "LXXXIV"),
    (85, "LXXXV"),
    (86, "LXXXVI"),
    (87, "LXXXVII"),
    (88, "LXXXVIII"),
    (89, "LXXXIX"),
    (90, "XC"),
    (91, "XCI"),
    (92, "XCII"),
    (93, "XCIII"),
    (94, "XCIV"),
    (95, "XCV"),
    (96, "XCVI"),
    (97, "XCVII"),
    (98, "XCVIII"),
    (99, "XCIX"),
    (100, "C"),
    (101, "CI"),
    (102, "CII"),
    (103, "CIII"),
    (104, "CIV"),
    (105, "CV"),
    (106, "CVI"),
    (107, "CVII"),
    (108, "CVIII"),
    (109, "CIX"),
    (110, "CX"),
    (111, "CXI"),
    (112, "CXII"),
    (113, "CXIII"),
    (114, "CXIV"),
    (115, "CXV"),
    (116, "CXVI"),
    (117, "CXVII"),
    (118, "CXVIII"),
    (119, "CXIX"),
    (120, "CXX"),
    (121, "CXXI"),
    (122, "CXXII"),
    (188, "CLXXXVIII"),
    (189, "CLXXXIX"),
    (198, "CXCVIII"),
    (199, "CXCIX"),
    (200, "CC"),
    (9, "IX"),
    (np.int16(10), "X"),
    (np.int32(11), "XI"),
    (np.int64(14), "XIV"),
]

exceptions_table = [
    (roman.to_roman, "X", TypeError),
    (roman.to_roman, -1, roman.OutOfRangeError),
    (roman.to_roman, 0, roman.OutOfRangeError),
    (roman.to_roman, 5000, roman.OutOfRangeError),
    (roman.from_roman, "asdfasd", roman.InvalidRomanNumeralError),
    (roman.from_roman, 1, TypeError),
    (roman.from_roman, "xi", roman.InvalidRomanNumeralError),
    (roman.is_roman_numeral, 1, TypeError),  # TODO: tbh I would just return False here?
    (roman.is_roman_numeral, ("I", "II"), TypeError),
]


@pytest.mark.parametrize("integer, roman_numeral", ints_and_roman_numerals)
def test_to_roman(integer, roman_numeral):
    """
    Test that `~plasmapy.utils.roman.to_roman` correctly converts
    integers to Roman numerals, and that the inverse is true as well.
    """
    run_test(func=roman.to_roman, args=integer, expected_outcome=roman_numeral)
    run_test(func=roman.from_roman, args=roman_numeral, expected_outcome=int(integer))


@pytest.mark.parametrize("function, argument, expected_exception", exceptions_table)
def test_to_roman_exceptions(function, argument, expected_exception):
    """
    Test that `~plasmapy.utils.roman` functions raise the correct
    exceptions when necessary.
    """
    run_test(func=function, args=argument, expected_outcome=expected_exception)


test_is_roman_numeral_table = [
    ("I", True),
    ("i", False),
    ("CLXXXVIII", True),
    ("khjfda", False),
    ("VIIII", False),
    ("IXX", False),
]


@pytest.mark.parametrize("argument, expected", test_is_roman_numeral_table)
def test_is_roman_numeral(argument, expected):
    run_test(func=roman.is_roman_numeral, args=argument, expected_outcome=expected)

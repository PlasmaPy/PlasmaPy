import unittest
from collections import defaultdict
from plasmapy.utils import compare_default_dicts


class TestDictCompare(unittest.TestCase):
    """
    Tests custom comparing of pythons defaultdicts.
    Note that default value of int = 0 and default of
    str = ''
    """

    def test_same_keys(self):
        msg = "Two equal defaultdicts should be same."
        a = defaultdict(int, a=1, b=1, c=1)
        b = defaultdict(int, a=1, b=1, c=1)
        result_a_first = compare_default_dicts(a, b)
        result_b_first = compare_default_dicts(b, a)
        self.assertEqual(result_a_first, True, msg)
        self.assertEqual(result_b_first, True, msg)

    def test_different_values_on_intersecting_keys_same_default(self):
        msg = "Same key, same default, but different value should fail."
        a = defaultdict(int, a=2, b=2, c=2)
        b = defaultdict(int, a=1, b=1, c=1)
        result_a_first = compare_default_dicts(a, b)
        result_b_first = compare_default_dicts(b, a)
        self.assertEqual(result_a_first, False, msg)
        self.assertEqual(result_b_first, False, msg)

    @unittest.skip("Awaiting specific requirements")
    def test_disjoint_values_match_default_of_other(self):
        msg = "???"
        a = defaultdict(str, a=0, b=0, c="")
        b = defaultdict(int, c="", d="", e="")
        result_a_first = compare_default_dicts(a, b)
        result_b_first = compare_default_dicts(b, a)
        self.assertEqual(result_a_first, True, msg)
        self.assertEqual(result_b_first, True, msg)

    def test_two_empty_with_same_default_value(self):
        msg = "Empty defaultdicts are equal."
        a = defaultdict(int)
        b = defaultdict(int)
        result_a_first = compare_default_dicts(a, b)
        result_b_first = compare_default_dicts(b, a)
        self.assertEqual(result_a_first, True, msg)
        self.assertEqual(result_b_first, True, msg)

    @unittest.skip("Awaiting specific requirements")
    def test_two_empty_with_different_default_value(self):
        msg = "???"
        a = defaultdict(int)
        b = defaultdict(str)
        result_a_first = compare_default_dicts(a, b)
        result_b_first = compare_default_dicts(b, a)
        self.assertEqual(result_a_first, True, msg)
        self.assertEqual(result_b_first, True, msg)

    def test_no_intersecting_keys_same_default_value(self):
        msg = "No intersecting keys but same default value."
        a = defaultdict(int, a=0, b=0, c=0)
        b = defaultdict(int, d=0, e=0, f=0)
        result_a_first = compare_default_dicts(a, b)
        result_b_first = compare_default_dicts(b, a)
        self.assertEqual(result_a_first, True, msg)
        self.assertEqual(result_b_first, True, msg)

    def test_no_intersecting_keys_different_default_value(self):
        msg = "No intersecting keys but different disjoint \
               values AND different default value are NOT equal."
        a = defaultdict(str, a='', b='', c='')
        b = defaultdict(int, d=0, e=0, f=0)
        result_a_first = compare_default_dicts(a, b)
        result_b_first = compare_default_dicts(b, a)
        self.assertEqual(result_a_first, False, msg)
        self.assertEqual(result_b_first, False, msg)

    def test_some_intersecting_same_default_different_disjoint_values(self):
        msg = "Keys unique to one should yield default value if \
               queried on other."
        a = defaultdict(int, a=10, b=10, c=5)
        b = defaultdict(int, c=5, d=20, e=20)
        result_a_first = compare_default_dicts(a, b)
        result_b_first = compare_default_dicts(b, a)
        self.assertEqual(result_a_first, False, msg)
        self.assertEqual(result_b_first, False, msg)

    def test_comparing_should_not_mutate_any_dict(self):
        msg = "Comparing should not change any dict."
        a = defaultdict(int, a='Hello', b=0.5, c=[])
        a_copy = a.copy()
        b = defaultdict(str, a=1001, b={}, c='World', d=[])
        b_copy = b.copy()

        # We're not interested in the actual dictionaries,
        # just whether comparing mutated any of them.
        _ = compare_default_dicts(a, b)
        self.assertEqual(a.items(), a_copy.items(), msg)
        self.assertEqual(b.items(), b_copy.items(), msg)

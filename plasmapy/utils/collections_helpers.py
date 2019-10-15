from collections import defaultdict


def compare_default_dicts(a: defaultdict, b: defaultdict) -> bool:
    """Compare two defaultdicts, return True if equal, else False.
    Does a benign or soft compare. If the defaultdicts COULD become
    equal, they are considered equal.

    * Does NOT change the memory imprint of any of the dictionaries.
    * Any overlapping keys, must have same value.
    * Keys unique to one, must have the default value of the other.
    * Order of input does NOT matter.

    Example:
    a = defaultdict(lambda: "", a=42, b=42, c="")
    b = defaultdict(lambda: 42, c="", d="")
    compare_defaultdicts(a, b) -> True

    Parameters
    ----------
    a : defaultdict
        Default dictionary from collections
    b : defaultdict:
        Default dictionary from collections

    Returns
    -------
    bool : True if equal, else False

    """
    a_keys = set(a)
    b_keys = set(b)
    a_unique_keys = (a_keys | b_keys) - b_keys
    b_unique_keys = (a_keys | b_keys) - a_keys

    # The intersecting keys must have the same value
    if not all(a[key] == b[key] for key in (a_keys & b_keys)):
        return False

    # Keys unique to one, must have default value of other.
    if not all(b.default_factory() == a[key] for key in a_unique_keys):
        return False
    if not all(a.default_factory() == b[key] for key in b_unique_keys):
        return False

    return True

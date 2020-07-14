import json
from plasmapy.particles.particle_class import ParticleJSONDecoder


def json_load_particle(
    fp,
    *,
    cls=ParticleJSONDecoder,
    object_hook=None,
    parse_float=None,
    parse_int=None,
    parse_constant=None,
    object_pairs_hook=None,
    **kw
):
    """
    Returns the appropriate object from input JSON file object
    For more information, refer to documentation:
    https://docs.python.org/3/library/json.html#json.load

    Parameters
    ----------
    fp: file object
        file object of the JSON file that contains the JSON string describing
        the particle
    object_hook:
        function called with result of decoded literal.
        Refer to link above for more information.
    parse_float:
        function applied when decoding floats.
        Refer to link above for more information.
    parse_int:
        function applied when decoding integers.
        Refer to link above for more information.
    parse_constant:
        function applied when decoding constants.
        Refer to link above for more information.
    object_pairs_hook:
        function called with result of decode with ordered list of pairs
        Refer to link above for more information.
    """
    return json.load(
        fp,
        cls=cls,
        object_hook=object_hook,
        parse_float=parse_float,
        parse_int=parse_int,
        parse_constant=parse_constant,
        object_pairs_hook=object_pairs_hook,
    )


def json_loads_particle(
    s,
    *,
    cls=ParticleJSONDecoder,
    object_hook=None,
    parse_float=None,
    parse_int=None,
    parse_constant=None,
    object_pairs_hook=None,
    **kw
):
    """
    Returns the appropriate object from input JSON string representation
    For more information refer to documentation:
    https://docs.python.org/3/library/json.html#json.loads

    Parameters
    ----------
    s: string
        JSON string literal to be deserialized into a particle
    object_hook:
        function called with result of decoded literal.
        Refer to link above for more information.
    parse_float:
        function applied when decoding floats.
        Refer to link above for more information.
    parse_int:
        function applied when decoding integers.
        Refer to link above for more information.
    parse_constant:
        function applied when decoding constants.
        Refer to link above for more information.
    object_pairs_hook:
        function called with result of decode with ordered list of pairs
        Refer to link above for more information.
    """
    return json.loads(
        s,
        cls=cls,
        object_hook=object_hook,
        parse_float=parse_float,
        parse_int=parse_int,
        parse_constant=parse_constant,
        object_pairs_hook=object_pairs_hook,
    )

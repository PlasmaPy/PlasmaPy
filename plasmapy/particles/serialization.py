__all__ = [
    "json_load_particle",
    "json_loads_particle",
    "ParticleJSONDecoder",
]

import json
from plasmapy.particles.exceptions import InvalidElementError
from plasmapy.particles.particle_class import (
    AbstractParticle,
    CustomParticle,
    DimensionlessParticle,
    Particle,
)


class ParticleJSONDecoder(json.JSONDecoder):
    """
    A custom `~json.JSONDecoder` class to deserialize JSON objects into
    PlasmaPy particle objects.

    Parameters
    ----------
    object_hook :
        If specified, will be called with the result of every JSON object
        decoded and its return value will be used in place of the given `dict`.
        This can be used to provide custom deserializations (e.g. to support
        JSON-RPC class hinting)
        (If not specified, then defaults to `particle_hook`.).
    **kwargs :
        Any keyword accepted by `~json.JSONDecoder`.
    """

    def __init__(self, *, object_hook=None, **kwargs):
        if object_hook is None:
            object_hook = self.particle_hook
        json.JSONDecoder.__init__(self, object_hook=object_hook, **kwargs)

    @staticmethod
    def particle_hook(json_dict):
        """
        An `object_hook` utilized by the `json` deserialization processes to decode
        json strings into a `plasmapy` particle class (`AbstractParticle`,
        `CustomParticle`, `DimensionlessParticle`, `Particle`).
        """
        particle_types = {
            "AbstractParticle": AbstractParticle,
            "CustomParticle": CustomParticle,
            "DimensionlessParticle": DimensionlessParticle,
            "Particle": Particle,
        }
        if "plasmapy_particle" in json_dict:
            try:
                pardict = json_dict["plasmapy_particle"]
                partype = pardict["type"]
                args = pardict["__init__"]["args"]
                kwargs = pardict["__init__"]["kwargs"]
                particle = particle_types[partype](*args, **kwargs)
                return particle
            except KeyError:
                raise InvalidElementError(
                    f"json file does not define a valid plasmapy particle"
                )
        else:
            return json_dict


def json_load_particle(fp, *, cls=ParticleJSONDecoder, **kwargs):
    """
    A convenient form of `json.load` to deserialize a JSON document into a
    PlasmaPy particle object. (Mirrors `json.load` with `cls` defaulting to
    `ParticleJSONDecoder`.).

    Parameters
    ----------
    fp : `file object <https://docs.python.org/3/glossary.html#term-file-object>`_
        A file object containing a JSON document.
    cls : `json.JSONDecoder` class
        A `~json.JSONDecoder` class. (Default `ParticleJSONDecoder`).
    **kwargs :
        Any keyword accepted by `json.load`.
    """
    return json.load(fp, cls=cls, **kwargs)


def json_loads_particle(s, *, cls=ParticleJSONDecoder, **kwargs):
    """
    A convenient form of `json.loads` to deserialize a JSON string into a
    PlasmaPy particle object. (Mirrors `json.loads` with `cls` defaulting to
    `ParticleJSONDecoder`.).

    Parameters
    ----------
    s : str
        A JSON string.
    cls : `json.JSONDecoder` class
        A `~json.JSONDecoder` class. (Default `ParticleJSONDecoder`).
    **kwargs :
        Any keyword accepted by `json.loads`.
    """
    return json.loads(s, cls=cls, **kwargs)

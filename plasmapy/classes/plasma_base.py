from collections import OrderedDict
from abc import ABC

# GenericPlasma subclass registry
PLASMA_CLASSES = OrderedDict()

class GenericPlasmaRegistrar(ABC):
    """
    Registration class for `~plasmapy.classes.GenericPlasma`.

    This class checks for the existance of a method named ``is_datasource_for``
    when a subclass of `GenericPlasma` is defined. If it exists it will add that
    class to the registry.
    """
    _registry = PLASMA_CLASSES

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if hasattr(cls, 'is_datasource_for'):
            cls._registry[cls] = cls.is_datasource_for


class GenericPlasma(GenericPlasmaRegistrar):
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

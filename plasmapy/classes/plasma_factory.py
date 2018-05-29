from plasmapy.classes.plasma_base import GenericPlasma

from plasmapy.utils.datatype_factory_base import BasicRegistrationFactory
from plasmapy.utils.datatype_factory_base import NoMatchError
from plasmapy.utils.datatype_factory_base import MultipleMatchError
from plasmapy.utils.datatype_factory_base import ValidationFunctionError

__all__ = [
    "PlasmaFactory",
    "Plasma",
]

class PlasmaFactory(BasicRegistrationFactory):
    """
    Plasma factory class. Used to create a variety of Plasma objects.
    Valid plasma structures are specified by registering them with the
    factory.
    """
    pass


Plasma = PlasmaFactory(default_widget_type=GenericPlasma,
                       registry=GenericPlasma._registry,
                       additional_validation_functions=['is_datasource_for'])

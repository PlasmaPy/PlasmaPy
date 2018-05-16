from plasmapy.classes.plasma_base import GenericPlasma, PLASMA_CLASSES

from plasmapy.utils.datatype_factory_base import BasicRegistrationFactory
from plasmapy.utils.datatype_factory_base import NoMatchError
from plasmapy.utils.datatype_factory_base import MultipleMatchError
from plasmapy.utils.datatype_factory_base import ValidationFunctionError


class PlasmaFactory(BasicRegistrationFactory):
    """
    Plasma factory class. Used to create a variety of Map objects. Valid plasma
    structures are specified by registering them with the factory.
    """
    pass


Plasma = PlasmaFactory(default_widget_type=GenericPlasma,
                       additional_validation_functions=['is_datasource_for'])
Plasma.registry = PLASMA_CLASSES

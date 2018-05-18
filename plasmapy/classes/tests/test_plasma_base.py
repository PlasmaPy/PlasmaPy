from plasmapy.classes.plasma_base import (
                        GenericPlasmaRegistrar,
                        GenericPlasma,
                        )

class NoDataSource(GenericPlasmaRegistrar):
    pass

class IsDataSource(GenericPlasmaRegistrar):
    @classmethod
    def is_datasource_for(cls, **kwargs):
        return True

class IsNotDataSource(GenericPlasmaRegistrar):
    @classmethod
    def is_datasource_for(cls, **kwargs):
        return False


class TestRegistrar:
    def test_no_data_source(self):
        """
        NoDataSource class should not be registered since it has
        no method named ``is_datasource_for``.
        """
        assert not GenericPlasmaRegistrar._registry.get(NoDataSource)

    def test_is_data_source(self):
        """
        IsDataSource class should be registered since it has a
        method named ``is_datasource_for`` and must return True.
        """
        assert GenericPlasmaRegistrar._registry.get(IsDataSource)
        assert GenericPlasmaRegistrar._registry[IsDataSource]()
        # Delete the class from _registry once test is done
        # to not interfere with plasma factory tests
        del GenericPlasmaRegistrar._registry[IsDataSource]

    def test_is_not_data_source(self):
        """
        IsNotDataSource class should be registered since it has a
        method named ``is_datasource_for`` but must return False.
        """
        assert GenericPlasmaRegistrar._registry.get(IsNotDataSource)
        assert not GenericPlasmaRegistrar._registry[IsNotDataSource]()
        del GenericPlasmaRegistrar._registry[IsNotDataSource]


def test_subclasses():
    assert issubclass(GenericPlasma, GenericPlasmaRegistrar)

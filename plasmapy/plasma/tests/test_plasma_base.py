from plasmapy.plasma import BasePlasma, GenericPlasma

# Get rid of any previously registered classes.
BasePlasma._registry = {}


class NoDataSource(BasePlasma):
    pass


class IsDataSource(BasePlasma):
    @classmethod
    def is_datasource_for(cls, **kwargs):
        return True


class IsNotDataSource(BasePlasma):
    @classmethod
    def is_datasource_for(cls, **kwargs):
        return False


class TestRegistrar:
    def test_no_data_source(self):
        """
        NoDataSource class should not be registered since it has
        no method named ``is_datasource_for``.
        """
        assert not BasePlasma._registry.get(NoDataSource)

    def test_is_data_source(self):
        """
        IsDataSource class should be registered since it has a
        method named ``is_datasource_for`` and must return True.
        """
        assert BasePlasma._registry.get(IsDataSource)
        assert BasePlasma._registry[IsDataSource]()
        # Delete the class from _registry once test is done
        # to not interfere with plasma factory tests
        del BasePlasma._registry[IsDataSource]

    def test_is_not_data_source(self):
        """
        IsNotDataSource class should be registered since it has a
        method named ``is_datasource_for`` but must return False.
        """
        assert BasePlasma._registry.get(IsNotDataSource)
        assert not BasePlasma._registry[IsNotDataSource]()
        del BasePlasma._registry[IsNotDataSource]


def test_subclasses():
    assert issubclass(GenericPlasma, BasePlasma)

.. _classes:

===============
PlasmaPy Plasma
===============

.. module:: plasmapy.classes
.. currentmodule:: plasmapy.classes

Overview
--------
One of the core classes in PlasmaPy is `~plasmapy.classes.Plasma`. In order
to make it easy to work with different plasma data in PlasmaPy, the
`~plasmapy.classes.Plasma` object provides a number of methods for
commonly existing plasmas in nature.

All Plasma objects are created using the Plasma factory
`~plasmapy.classes.Plasma`.

A number of plasma data structures are supported by subclassing this
base object. See :ref:`plasma-sources` to see a list of all of them.

Creating Plasma Objects
-----------------------
Plasma objects are constructed using the special factory class
`~plasmapy.classes.Plasma`: ::

    >>> x = plasmapy.classes.Plasma(T_e=T_e,
    ...                             n_e=n_e,
    ...                             Z=Z,
    ...                             particle=particle)  # doctest: +SKIP

The result of a call to `~plasmapy.classes.Plasma` will be either a
`~plasmapy.classes.GenericPlasma` object, or a subclass of
`~plasmapy.classes.GenericPlasma` which deals with a specific type of
data, e.g. `~plasmapy.classes.sources.PlasmaBlob` or
`~plasmapy.classes.sources.Plasma3D` (see :ref:`plasma-sources` to see
a list of all of them).

.. autoclass:: plasmapy.classes.plasma_factory.PlasmaFactory

Using Plasma Objects
--------------------

Once a Plasma object has been created using `~plasmapy.classes.Plasma`
it will be a instance or a subclass of the
`~plasmapy.classes.GenericPlasma` class. The documentation of
`~plasmapy.classes.GenericPlasma` lists the attributes and methods that
are available on all Plasma objects.

.. _plasma-classes:

Plasma Classes
--------------
Defined in `plasmapy.classes.sources` are a set of
`~plasmapy.classes.GenericPlasma` subclassesmwhich convert the keyword
arguments data to the standard `~plasmapy.classes.GenericPlasma`
interface. These subclasses also provide a method, which describes to
the `Plasma <plasmapy.classes.plasma_factory.PlasmaFactory>` factory
which the data match its plasma data structure.

.. automodapi:: plasmapy.classes
    :no-main-docstr:
    :no-heading:

.. _plasma-sources:

Plasma Subclasses
-----------------

.. toctree::
    sources

Writing a new Plasma subclass
--------------------------------

Any subclass of `~plasmapy.classes.GenericPlasma` which defines a method
named `~plasmapy.classes.GenericPlasma.is_datasource_for` will
automatically be registered with the
`Plasma <plasmapy.classes.plasma_factory.PlasmaFactory>` factory. The
``is_datasource_for`` method describes the form of the data for which
the `~plasmapy.classes.GenericPlasma` subclass is valid. For example,
it might check the number and types of keyword arguments. This makes it
straightforward to define your own `~plasmapy.classes.GenericPlasma`
subclass for a new data structure or a custom data source like simulated
data. These classes only have to be imported for this to work, as
demonstrated by the following example.

.. code-block:: python

    import plasmapy.classes
    import astropy.units as u

    class FuturePlasma(plasmapy.classes.GenericPlasma):
        def __init__(self, **kwargs):

            super(FuturePlasma, self).__init__(**kwargs)

        # Specify a classmethod that determines if the input data matches
        # this new subclass
        @classmethod
        def is_datasource_for(cls, **kwargs):
            """
            Determines if any of keyword arguments have a dimensionless value.
            """
            for _, value in kwargs.items():
                try:
                    if value.unit == u.dimensionless_unscaled:
                        return True
                except AttributeError:
                    pass

            return False


This class will now be available through the
`Plasma <plasmapy.classes.plasma_factory.PlasmaFactory>` factory as long
as this class has been defined, i.e. imported into the current session.

If you do not want to create a method named ``is_datasource_for`` you
can manually register your class and matching method using the following
method.

.. code-block:: python

    import plasmapy.classes

    plasmapy.classes.Plasma.register(FuturePlasma, FuturePlasma.some_matching_method)

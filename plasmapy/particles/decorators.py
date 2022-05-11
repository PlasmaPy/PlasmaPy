"""Decorators for `plasmapy.particles`."""

__all__ = ["ParticleValidator", "particle_input"]

import functools
import inspect
import wrapt

from numbers import Integral
from typing import (
    Any,
    Callable,
    Dict,
    List,
    NoReturn,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
)

from plasmapy.particles._factory import _physical_particle_factory
from plasmapy.particles.exceptions import (
    ChargeError,
    InvalidElementError,
    InvalidIonError,
    InvalidIsotopeError,
    ParticleError,
)
from plasmapy.particles.particle_class import Particle, ParticleLike
from plasmapy.particles.particle_collections import ParticleList

# Temporarily define ParticleListLike, pending #1528
ParticleListLike = Union[ParticleList, Sequence[ParticleLike]]


def _get_annotations(f: Callable):
    """
    Access the annotations of a callable `object`.

    Notes
    -----
    This function should be replaced with inspect.get_annotations when
    support for Python 3.9 is dropped.
    """
    return getattr(f, "__annotations__", None)


def _make_into_set_or_none(obj) -> Optional[Set]:
    """
    Return `None` is ``obj`` is `None`, and otherwise convert ``obj``
    into a `set`.

    If ``obj`` is a string, then a `set` containing only ``obj`` will
    be returned (i.e., ``obj`` will not be treated as iterable).
    """
    if obj is None:
        return obj
    return {obj} if isinstance(obj, str) else set(obj)


def _bind_arguments(
    func: Callable,
    args: Optional[Tuple] = None,
    kwargs: Optional[Dict[str, Any]] = None,
    instance=None,
):
    """
    Bind the arguments provided by ``args`` and ``kwargs`` to
    the corresponding parameters in the signature of the function
    or method being decorated.

    Parameters
    ----------
    func : callable
        The function or method to which to bind ``args`` and ``kwargs``.

    args : tuple
        Positional arguments.

    kwargs : dict
        Keyword arguments.

    instance
        If ``func`` is a class instance method, then ``instance`` should
        be the instance to which ``func`` belongs.

    Returns
    -------
    dict
        A dictionary with the parameters of the wrapped function
        as keys and the corresponding arguments as values,
        but removing ``self`` and ``cls``.
    """
    wrapped_signature = inspect.signature(func)

    # When decorating a function or staticmethod, instance will
    # be None. When decorating a class instance method, instance
    # will be the class instance, and will need to be bound to
    # the "self" parameter but later removed. For a class method,
    # it will be bound to the "cls" parameter instead.

    if instance is None:
        bound_arguments = wrapped_signature.bind(*args, **kwargs)
    else:
        bound_arguments = wrapped_signature.bind(instance, *args, **kwargs)

    bound_arguments.apply_defaults()
    arguments_to_be_processed = bound_arguments.arguments

    arguments_to_be_processed.pop("self", None)
    arguments_to_be_processed.pop("cls", None)

    return arguments_to_be_processed


_basic_allowed_annotations = (
    Particle,  # deprecated
    ParticleLike,
    ParticleListLike,
    Union[ParticleLike, ParticleListLike],
    (Particle, Particle),  # deprecated
)
_optional_allowed_annotations = tuple(
    Optional[annotation]
    for annotation in _basic_allowed_annotations
    if annotation != (Particle, Particle)  # temporary hack
)
_allowed_annotations = _basic_allowed_annotations + _optional_allowed_annotations


class ParticleValidator:
    """
    Processes arguments for |particle_input|.

    Parameters
    ----------
    wrapped : function
        The function or method to be decorated.

    require : `str`, `set`, `list`, or `tuple`, optional
        Categories that a particle must be in.  If a particle is not in
        all of these categories, then a
        `~plasmapy.particles.exceptions.ParticleError` will be raised.

    any_of : `str`, `set`, `list`, or `tuple`, optional
        Categories that a particle may be in.  If a particle is not in
        any of these categories, then a
        `~plasmapy.particles.exceptions.ParticleError` will be raised.

    exclude : `str`, `set`, `list`, or `tuple`, optional
        Categories that a particle cannot be in.  If a particle is in
        any of these categories, then a
        `~plasmapy.particles.exceptions.ParticleError` will be raised.

    allow_custom_particles : bool
        If `True`, allow |CustomParticle| instances to be passed through.
        Defaults to `True`.
    """

    def __init__(
        self,
        wrapped: Callable,
        *,
        require: Optional[Union[str, Set, List, Tuple]] = None,
        any_of: Optional[Union[str, Set, List, Tuple]] = None,
        exclude: Optional[Union[str, Set, List, Tuple]] = None,
        allow_custom_particles: bool = True,
    ):
        self._data = {}
        self.wrapped = wrapped
        self.require = require
        self.any_of = any_of
        self.exclude = exclude
        self.allow_custom_particles = allow_custom_particles

    @property
    def wrapped(self) -> Callable:
        """
        The function that is being decorated.

        Returns
        -------
        callable
        """
        return self._data["wrapped_function"]

    @wrapped.setter
    def wrapped(self, function: Callable):
        self._data["wrapped_function"] = function
        self._data["annotations"] = _get_annotations(function)
        self._data["parameters_to_process"] = self.find_parameters_to_process()

    def find_parameters_to_process(self) -> List[str]:
        """
        Identify the parameters that have annotations that indicate that
        the

        Returns
        -------
        `list` of `str`
        """
        return [
            parameter
            for parameter, annotation in self.annotations.items()
            if annotation in _allowed_annotations and parameter != "return"
        ]

    @property
    def annotations(self) -> Dict[str, Any]:
        """
        The annotations of the decorated function.

        Returns
        -------
        dict
        """
        return self._data.get("annotations", None)

    @property
    def require(self) -> Optional[Set]:
        """
        Categories that the particle must belong to.

        Returns
        -------
        `set` of `str` or `None`
        """
        return self._data["require"]

    @require.setter
    def require(self, require_: Optional[Union[str, Set, List, Tuple]]):
        self._data["require"] = _make_into_set_or_none(require_)

    @property
    def any_of(self) -> Optional[Set]:
        """
        Categories of which the particle must belong to at least one.

        Returns
        -------
        `set` of `str` or `None`
        """
        return self._data["any_of"]

    @any_of.setter
    def any_of(self, any_of_: Optional[Union[str, Set, List, Tuple]]):
        self._data["any_of"] = _make_into_set_or_none(any_of_)

    @property
    def exclude(self) -> Optional[Set]:
        """
        Categories that the particle cannot belong to.

        Returns
        -------
        `set` of `str` or `None`
        """
        return self._data["exclude"]

    @exclude.setter
    def exclude(self, exclude_):
        self._data["exclude"] = _make_into_set_or_none(exclude_)

    @property
    def allow_custom_particles(self) -> bool:
        """
        If `True`,  then the decorated argument may be or include
        |CustomParticle| instances. Defaults to `True`.

        Returns
        -------
        bool
        """
        return self._data["allow_custom_particles"]

    @allow_custom_particles.setter
    def allow_custom_particles(self, allow_custom_particles_: bool):
        self._data["allow_custom_particles"] = allow_custom_particles_

    @property
    def parameters_to_process(self) -> List[str]:
        """
        The parameters of
        `~plasmapy.particles.decorators2.ParticleValidator.wrapped`
        that have annotations to be processed by |particle_input|.

        Returns
        -------
        `list` of `str`
        """
        return self._data["parameters_to_process"]

    def _verify_charge_categorization(self, particle) -> NoReturn:
        """
        Raise an exception if the particle is required to have charge
        information and does not, or if the particle is required to be
        charged and is not.

        Raises
        ------
        ~plasmapy.particles.exceptions.ChargeError
            If the particle is required to have charge information and
            does not, or if the particle is required to be charged and
            is either uncharged or lacks charge information.
        """
        must_be_charged = self.require is not None and "charged" in self.require
        must_have_charge_info = self.any_of == {"charged", "uncharged"}

        uncharged = particle._attributes["charge number"] == 0
        lacks_charge_info = particle._attributes["charge number"] is None

        if must_be_charged and (uncharged or must_have_charge_info):
            raise ChargeError(f"A charged particle is required for {self.wrapped}.")

        if must_have_charge_info and lacks_charge_info:
            raise ChargeError(f"Charge information is required for {self.wrapped}.")

    @staticmethod
    def _category_errmsg(particle, require, exclude, any_of, funcname) -> str:
        """
        Return an appropriate error message for when a particle does not
        meet the required categorical specifications.

        Returns
        -------
        str
        """
        category_errmsg = (
            f"The particle {particle} does not meet the required "
            f"classification criteria to be a valid input to {funcname}. "
        )

        errmsg_table = [
            (require, "must belong to all"),
            (any_of, "must belong to any"),
            (exclude, "cannot belong to any"),
        ]

        for condition, phrase in errmsg_table:
            if condition:
                category_errmsg += (
                    f"The particle {phrase} of the following categories: {condition}. "
                )

        return category_errmsg

    def _verify_particle_categorization(self, particle) -> NoReturn:
        """
        Verify that the particle meets the categorization criteria.

        Raises
        ------
        ~plasmapy.particles.exceptions.ParticleError
            If the particle does not meet the categorization criteria.

        See Also
        --------
        ~plasmapy.particles.particle_class.Particle.is_category
        """
        if not particle.is_category(
            require=self.require,
            any_of=self.any_of,
            exclude=self.exclude,
        ):
            errmsg = self._category_errmsg(
                particle,
                self.require,
                self.exclude,
                self.any_of,
                self.wrapped.__name__,
            )
            raise ParticleError(errmsg)

    def _verify_particle_name_criteria(self, parameter, particle):
        """
        Check that parameters with special names meet the expected
        criteria.
        """

        category_table = (
            ("element", particle.element, InvalidElementError),
            ("isotope", particle.isotope, InvalidIsotopeError),
            ("ion", particle.ionic_symbol, InvalidIonError),
        )

        for category_name, category_symbol, CategoryError in category_table:
            if parameter == category_name and not category_symbol:
                raise CategoryError(
                    f"The argument {parameter} = {parameter!r} to "
                    f"{self.wrapped.__name__} does not correspond to a "
                    f"valid {parameter}."
                )

    def process_argument(
        self,
        parameter: str,
        argument: Any,
        Z: Optional[Integral],
        mass_numb: Optional[Integral],
    ) -> Any:
        """
        Process an argument that has an appropriate annotation.

        Parameters
        ----------
        parameter : str
            The name of the :term:`parameter` that was decorated.

        argument : object
            The value of the :term:`argument` associated with ``parameter``.

        Z : integer
            The charge number of an ion or neutral particle.

        mass_numb : integer
            The mass number of an isotope.

        Raises
        ------
        CategoryError
        """
        annotation = self.annotations.get(parameter, None)

        if annotation not in _allowed_annotations:
            return argument

        if annotation in _optional_allowed_annotations and argument is None:
            return argument

        # This does not yet include cases like Optional[ParticleList],
        # Union[ParticleList, ParticleLike], etc. and thus needs updating.

        if annotation == (Particle, Particle):  # deprecated
            if not hasattr(argument, "__len__") or len(argument) != 2:
                raise TypeError(f"The length of {argument} must be 2.")
            return Particle(argument[0]), Particle(argument[1])

        if annotation in _basic_allowed_annotations and argument is None:
            raise TypeError(f"{parameter} may not be None.")

        particle = _physical_particle_factory(argument, Z=Z, mass_numb=mass_numb)

        self._verify_charge_categorization(particle)
        self._verify_particle_categorization(particle)
        self._verify_particle_name_criteria(parameter, particle)

        return particle

    parameters_to_skip = ("Z", "mass_numb")

    def _perform_pre_validations(self, Z, mass_numb):

        if not self.parameters_to_process:
            raise ParticleError(
                "No parameters have an annotation that will invoke particle_input."
            )

        Z_or_mass_numb = Z is not None or mass_numb is not None
        multiple_annotated_parameters = len(self.parameters_to_process) > 1

        if Z is not None and not isinstance(Z, Integral):
            raise TypeError("Z must be an integer.")

        if mass_numb is not None and not isinstance(mass_numb, Integral):
            raise TypeError("mass_numb must be an integer.")

        if Z_or_mass_numb and multiple_annotated_parameters:
            raise ParticleError(
                "The arguments Z and mass_numb are not allowed when more "
                "than one argument or keyword is annotated with Particle "
                "in functions decorated with particle_input."
            )

    def process_arguments(
        self, args: tuple, kwargs: Dict[str, Any], instance=None
    ) -> Dict[str, Any]:

        arguments = _bind_arguments(self.wrapped, args, kwargs, instance)

        Z = arguments.pop("Z", None)
        mass_numb = arguments.pop("mass_numb", None)

        self._perform_pre_validations(Z, mass_numb)

        return {
            parameter: self.process_argument(parameter, argument, Z, mass_numb)
            for parameter, argument in arguments.items()
        }


# TODO: add ionic_level naming scheme


def particle_input(
    wrapped_function=None,
    require: Union[str, Set, List, Tuple] = None,
    any_of: Union[str, Set, List, Tuple] = None,
    exclude: Union[str, Set, List, Tuple] = None,
    allow_custom_particles: bool = True,
):
    """
    Convert appropriately annotated arguments to particle objects.

    Take |particle-like| or |particle-list-like| arguments that are
    appropriately annotated (for example, with |ParticleLike| or
    |ParticleListLike|), and convert them to a |Particle|,
    |CustomParticle|, or |ParticleList|.

    If the particle(s) do not satisfy provided categorization criteria,
    then raise an appropriate exception. If the annotated parameter is
    named ``element``, ``isotope``, ``ion``, or ``ionic_level``, then
    the corresponding particle must be consistent with the name.

    The conversion of objects to a |Particle|, |CustomParticle|, or
    |ParticleList| is done in
    `~plasmapy.particles.decorators2.ParticleValidator`.

    Parameters
    ----------
    wrapped_function : `callable`
        The function or method to be decorated.

    require : `str`, `set`, `list`, or `tuple`, keyword-only, optional
        Categories that a particle must be in.  If a particle is not in
        all of these categories, then a |ParticleError| will be raised.

    any_of : `str`, `set`, `list`, or `tuple`, keyword-only, optional
        Categories that a particle may be in.  If a particle is not in
        any of these categories, then a |ParticleError| will be raised.

    exclude : `str`, `set`, `list`, or `tuple`, keyword-only, optional
        Categories that a particle cannot be in.  If a particle is in
        any of these categories, then a |ParticleError| will be raised.

    allow_custom_particles : bool, keyword-only, optional
        ...

    Notes
    -----
    If the annotated argument is named ``element``, ``isotope``, or ``ion``,
    then the decorator will raise an
    `~plasmapy.particles.exceptions.InvalidElementError`,
    `~plasmapy.particles.exceptions.InvalidIsotopeError`, or
    `~plasmapy.particles.exceptions.InvalidIonError` if the particle
    does not correspond to an element, isotope, or ion, respectively.

    If exactly one argument is annotated with
    `~plasmapy.particles.particle_class.ParticleLike`, then the keywords
    ``Z`` and ``mass_numb`` may be used to specify the charge number
    and/or mass number of an ion or isotope.  However, the decorated
    function must allow ``Z`` and/or ``mass_numb`` as keywords in order
    to enable this functionality.

    Raises
    ------
    `TypeError`
        If the annotated argument is not a `str`, `int`, `tuple`, `list`
        or `~plasmapy.particles.particle_class.Particle`; or if ``Z`` or
        ``mass_numb`` is not an `int`.

    `ValueError`
        If the number of input elements in a collection do not match the
        number of expected elements.

    `~plasmapy.particles.exceptions.InvalidParticleError`
        If the annotated argument does not correspond to a valid
        particle.

    `~plasmapy.particles.exceptions.InvalidElementError`
        If an annotated argument is named ``element``, and the input
        does not correspond to an element, isotope, or ion.

    `~plasmapy.particles.exceptions.InvalidIsotopeError`
        If an annotated argument is named ``isotope``, and the input
        does not correspond to an isotope or an ion of an isotope.

    `~plasmapy.particles.exceptions.InvalidIonError`
        If an annotated argument is named ``ion``, and the input does
        not correspond to an ion.

    `~plasmapy.particles.exceptions.ChargeError`
        If ``'charged'`` is in the ``require`` argument and the particle
        is not explicitly charged, or if ``any_of = {'charged',
        'uncharged'}`` and the particle does not have charge information
        associated with it.

    `~plasmapy.particles.exceptions.ParticleError`
        If an annotated argument does not meet the criteria set by the
        categories in the ``require``, ``any_of``, and ``exclude``
        keywords; if more than one argument is annotated and ``Z`` or
        ``mass_numb`` are used as arguments; or if none of the arguments
        have been annotated with `~plasmapy.particles.particle_class.Particle`.

    Examples
    --------
    The following simple decorated function returns the
    `~plasmapy.particles.particle_class.Particle` object created from
    the function's sole argument: ← needs updating

    .. code-block:: python

        from plasmapy.particles import particle_input, ParticleLike

        @particle_input
        def decorated_function(particle: ParticleLike):
            return particle

    This decorator may be used for methods in instances of classes, as
    in the following example:

    .. code-block:: python

        from plasmapy.particles import particle_input, ParticleLike

        class SampleClass:
            @particle_input
            def decorated_method(self, particle: ParticleLike):
                return particle

        sample_instance = SampleClass()
        sample_instance.decorated_method('Fe')

    Some functions may be intended to be used with only certain categories
    of particles.  The ``require``, ``any_of``, and ``exclude`` keyword
    arguments enable this functionality.

    .. code-block:: python

        from plasmapy.particles import particle_input, ParticleLike

        @particle_input(
            require={'matter'},
            any_of={'charged', 'uncharged},
            exclude={'neutrino', 'antineutrino'},
        )
        def selective_function(particle: ParticleLike):
            return particle

    * Discuss annotated class methods!  Include that ``@particle_input``
      should be the inner decorator and ``@classmethod`` should be the
      outer decorator in order for this to work correctly.
    """

    # The following pattern comes from the docs for wrapt, and requires
    # that the arguments to the decorator are keyword-only.
    if wrapped_function is None:
        return functools.partial(
            particle_input,
            require=require,
            any_of=any_of,
            exclude=exclude,
            allow_custom_particles=allow_custom_particles,
        )

    particle_validator = ParticleValidator(
        wrapped=wrapped_function,
        require=require,
        any_of=any_of,
        exclude=exclude,
        allow_custom_particles=allow_custom_particles,
    )

    @wrapt.decorator
    def wrapper(wrapped: Callable, instance, args: Tuple, kwargs: Dict[str, Any]):
        new_kwargs = particle_validator.process_arguments(args, kwargs, instance)
        return wrapped(**new_kwargs)

    return wrapper(wrapped_function)

"""Decorators for `plasmapy.particles`."""

__all__ = ["particle_input"]

import functools
import inspect
import wrapt

from numbers import Integral
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    NoReturn,
    Optional,
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
    InvalidParticleError,
    ParticleError,
)
from plasmapy.particles.particle_class import CustomParticle, Particle, ParticleLike
from plasmapy.particles.particle_collections import ParticleList, ParticleListLike

_basic_particle_input_annotations = (
    Particle,  # deprecated
    ParticleLike,
    ParticleListLike,
    Union[ParticleLike, ParticleListLike],
    (Particle, Particle),  # deprecated
)
_optional_particle_input_annotations = tuple(
    Optional[annotation]
    for annotation in _basic_particle_input_annotations
    if annotation != (Particle, Particle)  # temporary hack
)
_particle_input_annotations = (
    _basic_particle_input_annotations + _optional_particle_input_annotations
)


def _get_annotations(callable_: Callable):
    """
    Access the annotations of a callable.

    .. note::

       For Python 3.10+, this should be replaced with
       `inspect.get_annotations`.
    """
    # Python 3.10: Replace this with inspect.get_annotations
    return getattr(callable_, "__annotations__", None)


def _make_into_set_or_none(obj) -> Optional[Set]:
    """
    Return `None` if ``obj`` is `None`, and otherwise convert ``obj``
    into a `set`.

    If ``obj`` is a string, then a `set` containing only ``obj`` will
    be returned (i.e., ``obj`` will not be treated as iterable).
    """
    if obj is None:
        return obj
    return {obj} if isinstance(obj, str) else set(obj)


def _bind_arguments(
    callable_: Callable,
    args: Optional[Tuple] = None,
    kwargs: Optional[Dict[str, Any]] = None,
    instance=None,
) -> Dict:
    """
    Bind the arguments provided by ``args`` and ``kwargs`` to
    the corresponding parameters in the signature of the callable_
    or method being decorated.

    Parameters
    ----------
    callable_ : callable
        The function or method to which to bind ``args`` and ``kwargs``.

    args : tuple, optional
        Positional arguments.

    kwargs : `dict` of `str` to `object`, optional
        Keyword arguments.

    instance
        If ``callable_`` is a class instance method, then ``instance``
        should be the instance to which ``callable_`` belongs.

    Returns
    -------
    dict
        A dictionary with the parameters of ``callable_`` as keys and
        the corresponding arguments as values, but removing ``self`` and
        ``cls``.
    """
    wrapped_signature = inspect.signature(callable_)

    # When decorating a callable_ or staticmethod, instance will
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


class _ParticleInput:
    """
    Processes arguments for |particle_input|.

    Parameters
    ----------
    callable_ : callable
        The callable_ or method to be decorated.

    require : `str`, `set`, `list`, or `tuple`, optional, |keyword-only|
        Categories that a particle must be in.  If a particle is not in
        all of these categories, then a |ParticleError| will be raised.

    any_of : `str`, `set`, `list`, or `tuple`, optional, |keyword-only|
        Categories that a particle may be in.  If a particle is not in
        any of these categories, then a |ParticleError| will be raised.

    exclude : `str`, `set`, `list`, or `tuple`, optional, |keyword-only|
        Categories that a particle cannot be in.  If a particle is in
        any of these categories, then a |ParticleError| will be raised.

    allow_custom_particles : bool, optional, |keyword-only|, default: `True`
        If `True`, allow |CustomParticle| instances to be passed through.

    allow_particle_lists : bool, optional, |keyword-only|, default: `True`
        If `True`, allow |ParticleList| instances to be passed through.
    """

    def __init__(
        self,
        callable_: Callable,
        *,
        require: Optional[Union[str, Set, List, Tuple]] = None,
        any_of: Optional[Union[str, Set, List, Tuple]] = None,
        exclude: Optional[Union[str, Set, List, Tuple]] = None,
        allow_custom_particles: bool = True,
        allow_particle_lists: bool = True,
    ):
        self._data = {}
        self.callable_ = callable_
        self.require = require
        self.any_of = any_of
        self.exclude = exclude
        self.allow_custom_particles = allow_custom_particles
        self.allow_particle_lists = allow_particle_lists

    @property
    def callable_(self) -> Callable:
        """
        The callable that is being decorated.

        Returns
        -------
        callable
        """
        return self._data["callable"]

    @callable_.setter
    def callable_(self, callable_: Callable):
        self._data["callable"] = callable_
        self._data["annotations"] = _get_annotations(callable_)
        self._data["parameters_to_process"] = self.find_parameters_to_process()

    def find_parameters_to_process(self) -> List[str]:
        """
        Identify the parameters that have annotations to indicate that
        they should be processed.

        Returns
        -------
        `list` of `str`
        """
        return [
            parameter
            for parameter, annotation in self.annotations.items()
            if annotation in _particle_input_annotations and parameter != "return"
        ]

    @property
    def annotations(self) -> Dict[str, Any]:
        """
        The annotations of the decorated callable_.

        Returns
        -------
        `dict` of `str` to `object`
        """
        return self._data.get("annotations", None)

    @property
    def require(self) -> Optional[Set]:
        """
        Categories that the particle must belong to.

        Returns
        -------
        `set` of `str`, or `None`
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
        `set` of `str`, or `None`
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
        `set` of `str`, or `None`
        """
        return self._data["exclude"]

    @exclude.setter
    def exclude(self, exclude_):
        self._data["exclude"] = _make_into_set_or_none(exclude_)

    @property
    def allow_custom_particles(self) -> bool:
        """
        If `True`, then the decorated argument may be or include
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
    def allow_particle_lists(self) -> bool:
        """
        If `True`, then the decorated argument may be a |ParticleList|.
        Defaults to `True`.

        Returns
        -------
        bool
        """
        return self._data["allow_particle_lists"]

    @allow_particle_lists.setter
    def allow_particle_lists(self, allow_particle_lists_: bool):
        self._data["allow_particle_lists"] = allow_particle_lists_

    @property
    def parameters_to_process(self) -> List[str]:
        """
        The parameters of
        `~plasmapy.particles.decorators._ParticleInput.callable_` that have
        annotations to be processed by |particle_input|.

        Returns
        -------
        `list` of `str`
        """
        return self._data["parameters_to_process"]

    def verify_charge_categorization(self, particle) -> NoReturn:
        """
        Raise an exception if the particle does not meet charge
        categorization criteria.

        Raises
        ------
        |ChargeError|
            If the particle is required to have charge information and
            does not, or if the particle is required to be charged and
            is either uncharged or lacks charge information.
        """
        must_be_charged = self.require is not None and "charged" in self.require
        must_have_charge_info = self.any_of == {"charged", "uncharged"}

        uncharged = particle.is_category("uncharged")
        lacks_charge_info = particle.is_category(exclude={"charged", "uncharged"})

        if isinstance(uncharged, Iterable):
            uncharged = any(uncharged)
            lacks_charge_info = any(lacks_charge_info)

        if must_be_charged and (uncharged or must_have_charge_info):
            raise ChargeError(f"{self.callable_} can only accept charged particles.")

        if must_have_charge_info and lacks_charge_info:
            raise ChargeError(
                f"{self.callable_} can only accept particles which have "
                f"explicit charge information."
            )

    @staticmethod
    def category_errmsg(particle, require, exclude, any_of, callable_name) -> str:
        """
        Return an error message for when a particle does not meet
        categorization criteria.

        Returns
        -------
        str
        """
        category_errmsg = (
            f"The particle {particle} does not meet the required "
            f"classification criteria to be a valid input to {callable_name}. "
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

    def verify_particle_categorization(self, particle) -> NoReturn:
        """
        Verify that the particle meets the categorization criteria.

        Raises
        ------
        |ParticleError|
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
            errmsg = self.category_errmsg(
                particle,
                self.require,
                self.exclude,
                self.any_of,
                self.callable_.__name__,
            )
            raise ParticleError(errmsg)

    def verify_particle_name_criteria(self, parameter, particle):
        """
        Check that parameters with special names meet the expected
        categorization criteria.
        """

        name_categorization_exception = [
            ("element", {"require": "element"}, InvalidElementError),
            ("isotope", {"require": "isotope"}, InvalidIsotopeError),
            # TODO: In the future, "ion" should be changed to allow only
            # ions and not neutral atoms with an explicit charge.
            (
                "ion",
                {"require": "element", "any_of": {"charged", "uncharged"}},
                InvalidIonError,
            ),
            (
                "ionic_level",
                {"require": "element", "any_of": {"charged", "uncharged"}},
                InvalidIonError,
            ),
        ]

        for name, categorization, exception in name_categorization_exception:
            if parameter != name or particle is None:
                continue

            meets_name_criteria = particle.is_category(**categorization)

            if isinstance(particle, Iterable) and not isinstance(particle, str):
                meets_name_criteria = all(meets_name_criteria)

            if not meets_name_criteria:
                raise exception(
                    f"The argument {parameter} = {particle!r} to "
                    f"{self.callable_.__name__} does not correspond to a "
                    f"valid {parameter}."
                )

    def verify_allowed_types(self, particle):
        """
        Verify that the particle object contains only the allowed types
        of particles.
        """
        if not self.allow_custom_particles and isinstance(particle, CustomParticle):
            raise InvalidParticleError(
                f"{self.callable_.__name__} does not accept CustomParticle "
                f"or CustomParticle-like inputs."
            )

        if not self.allow_particle_lists and isinstance(particle, ParticleList):
            raise InvalidParticleError(
                f"{self.callable_.__name__} does not accept ParticleList "
                "or particle-list-like inputs."
            )

        if not self.allow_custom_particles and isinstance(particle, ParticleList):
            if any(particle.is_category("custom")):
                raise InvalidParticleError(
                    f"{self.callable_.__name__} does not accept CustomParticle "
                    f"or CustomParticle-like inputs."
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

        If the annotation is not one covered by |particle_input|, then
        this method will return ``argument`` without alteration.

        Otherwise, if the annotation is ``Optional[...]`` and
        ``argument`` is `None`, then this method will return `None`.

        Otherwise, this method will use ``argument`` (and possibly ``Z``
        and/or ``mass_numb``) to create a |Particle|, |CustomParticle|,
        or |ParticleList|. If the resulting particle object does not
        meet the specified criteria (if provided), then this method will
        raise an exception. Otherwise, the resulting particle object
        will be returned. See the docstring for |particle_input| for
        more details.

        Parameters
        ----------
        parameter : str
            The name of the |parameter| that was decorated.

        argument : object
            The value of the |argument| associated with ``parameter``.

        Z : integer, optional
            The |charge number| of an ion or neutral particle.

        mass_numb : integer, optional
            The mass number of an isotope.

        Returns
        -------
        object
            This method will return a |Particle|, |CustomParticle|,
            |ParticleList|, or `None` if the parameter has an annotation
            as described in the docstring for |particle_input|. For all
            other annotations, this method will return ``argument``
            without alteration.
        """
        annotation = self.annotations.get(parameter, None)

        if annotation not in _particle_input_annotations:
            return argument

        if annotation in _optional_particle_input_annotations and argument is None:
            return argument

        # This does not yet include cases like Optional[ParticleList],
        # Union[ParticleList, ParticleLike], etc. and thus needs updating.

        if annotation == (Particle, Particle):  # deprecated
            if not hasattr(argument, "__len__") or len(argument) != 2:
                raise ValueError(f"The length of {argument} must be 2.")
            return Particle(argument[0]), Particle(argument[1])

        if annotation in _basic_particle_input_annotations and argument is None:
            raise TypeError(f"{parameter} may not be None.")

        particle = _physical_particle_factory(argument, Z=Z, mass_numb=mass_numb)

        self.verify_charge_categorization(particle)
        self.verify_particle_categorization(particle)
        self.verify_particle_name_criteria(parameter, particle)
        self.verify_allowed_types(particle)

        return particle

    parameters_to_skip = ("Z", "mass_numb")

    def perform_pre_validations(self, Z, mass_numb):
        """
        Check that there are annotated parameters, that ``Z`` and
        ``mass_numb`` are integers, and that ``Z`` and ``mass_numb`` are
        not parameters when more than one parameter is annotated.
        """

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
                "in callables decorated with particle_input."
            )

    def process_arguments(
        self, args: tuple, kwargs: Dict[str, Any], instance=None
    ) -> Dict[str, Any]:
        """
        Process the arguments passed to the callable_ callable.

        Parameters
        ----------
        args : tuple
            Positional arguments passed to the callable_ callable.

        kwargs : `dict` of `str` to `object`
            Keyword arguments provided to the callable_ callable.

        instance : `object`, optional
            If the callable_ callable is a class instance method, then
            ``instance`` should be the class instance to which ``func``
            belongs.
        """

        arguments = _bind_arguments(self.callable_, args, kwargs, instance)

        Z = arguments.pop("Z", None)
        mass_numb = arguments.pop("mass_numb", None)

        self.perform_pre_validations(Z, mass_numb)

        return {
            parameter: self.process_argument(parameter, argument, Z, mass_numb)
            for parameter, argument in arguments.items()
        }


def particle_input(
    callable_: Optional[Callable] = None,
    *,
    require: Union[str, Set, List, Tuple] = None,
    any_of: Union[str, Set, List, Tuple] = None,
    exclude: Union[str, Set, List, Tuple] = None,
    allow_custom_particles: bool = True,
    allow_particle_lists: bool = True,
) -> Callable:
    """
    Convert |particle-like| |arguments| into particle objects.

    When a callable is |decorated| with |particle_input|,
    |particle-like| |arguments| that are appropriately |annotated|
    (i.e., with |ParticleLike| or |ParticleListLike|) will be converted
    into a |Particle|, |CustomParticle|, or |ParticleList|.

    The parameters ``Z`` and ``mass_numb`` may be used to specify the
    |charge number| of an ion and mass number of an isotope,
    respectively, as long as ``Z`` and/or ``mass_numb`` are |parameters|
    accepted by the callable and only one |parameter| is appropriately
    |annotated|.

    If the |annotation| is created using `~typing.Optional` (e.g.,
    ``Optional[ParticleLike]``), then `None` can be provided to
    ``callable_``.

    If the particle representation does not satisfy any categorization
    criteria that have been provided, then |particle_input| will raise
    an exception.

    If the |annotated| |parameter| is named ``element``, ``isotope``,
    ``ion``, or ``ionic_level``, then |particle_input| will raise an
    exception if the |argument| provided to the callable is not
    consistent with the |parameter|.

    .. note::

       |Annotated| |parameters| named ``ion`` and ``ionic_level`` accept
       neutral atoms as long as the charge number is explicitly defined.
       In the future, this functionality may change so that |parameters|
       named ``ion`` require a nonzero |charge number|.

    .. note::

       When both |particle_input| and |validate_quantities| are used to
       decorate a :term:`function`, they may be used in either order.
       When using both |particle_input| and |validate_quantities| to
       decorate a :term:`method`, |particle_input| should be the outer
       decorator and |validate_quantities| should be the inner
       decorator.

    .. note::

       When decorating a class method with |particle_input|,
       `classmethod` should be the outer decorator and |particle_input|
       should be the inner decorator, and the first argument
       (representing the class) must be named ``cls``. However,
       stacking `classmethod` and |particle_input| requires Python 3.9+.

    Parameters
    ----------
    callable_ : callable, optional
        The function or method to be decorated.

    require : `str`, `set`, `list`, or `tuple`, |keyword-only|, optional
        Categories that each particle are required to be in.

    any_of : `str`, `set`, `list`, or `tuple`, |keyword-only|, optional
        Categories of which each particle must belong to at least one.

    exclude : `str`, `set`, `list`, or `tuple`, |keyword-only|, optional
        Categories that each particle cannot be in.

    allow_custom_particles : bool, |keyword-only|, optional, default: `True`
        If `True`, allow |CustomParticle| instances to be passed through.

    allow_particle_lists : bool, |keyword-only|, optional, default: `True`
        If `True`, allow |ParticleList| instances to be passed through.

    Returns
    -------
    callable

    Raises
    ------
    `TypeError`
        If the annotated argument is not |particle-like|; or if ``Z`` or
        ``mass_numb`` is not an integer.

    |InvalidParticleError|
        If the annotated argument does not correspond to a valid
        particle, ``allow_custom_particles`` is `False` and the argument
        corresponds to a |CustomParticle|, or ``allow_particle_lists``
        is `False` and the argument corresponds to a |ParticleList|.

    |InvalidParticleError|
        If the decorated argument has charge and/or mass number
        information, and ``Z`` and/or ``mass_numb`` contain
        contradictory information.

    |InvalidElementError|
        If an annotated argument is named ``element``, and the input
        does not correspond to an element, isotope, or ion.

    |InvalidIsotopeError|
        If an annotated argument is named ``isotope``, and the input
        does not correspond to an isotope or an ion of an isotope.

    |InvalidIonError|
        If an annotated argument is named ``ion``, and the input does
        not correspond to an ion.

    |ChargeError|
        If ``"charged"`` is in the ``require`` argument and the particle
        is not explicitly charged, or if ``any_of = {"charged",
        "uncharged"}`` and the particle does not have charge information
        associated with it.

    |ParticleError|
        If the returned particle(s) do not meet the categorization
        criteria specified through ``require``, ``any_of``, or
        ``exclude``; or if none of the |parameters| of ``callable_``
        have been appropriately annotated.

    |UnitConversionError|
        If the annotated argument is a |Quantity|, but does not have a
        physical type of mass or charge.

    Warns
    -----
    |ParticleWarning|
        If decorated argument has charge and/or mass number information,
        and ``Z`` and/or ``mass_numb`` contain redundant information.

    See Also
    --------
    ~plasmapy.particles.particle_class.CustomParticle
    ~plasmapy.particles.particle_class.Particle
    ~plasmapy.particles.particle_collections.ParticleList
    ~plasmapy.utils.decorators.validators.validate_quantities

    Examples
    --------
    The |particle_input| decorator takes appropriately annotated
    |particle-like| or |particle-list-like| arguments that are provided
    to ``callable_``, and converts them into a |Particle|,
    |CustomParticle|, or |ParticleList| object.

    The following decorated function accepts a |particle-like| input and
    returns the corresponding particle object.

    >>> from plasmapy.particles import particle_input, ParticleLike
    >>> import astropy.units as u

    >>> @particle_input
    ... def get_particle(particle: ParticleLike):
    ...     return particle

    >>> get_particle("p+")
    Particle("p+")
    >>> get_particle(["p+", "e-"])
    ParticleList(['p+', 'e-'])
    >>> get_particle(1e-26 * u.kg)
    CustomParticle(mass=1e-26 kg, charge=nan C)

    If the annotation is constructed using `typing.Optional`, then the
    decorated callable will allow `None` to pass.

    >>> from typing import Optional
    >>> @particle_input
    ... def get_particle_or_none(particle: Optional[ParticleLike]):
    ...     return particle
    >>> get_particle_or_none("p+")
    Particle("p+")
    >>> print(get_particle_or_none(None))
    None

    If the decorated callable has parameters named ``Z`` and/or
    ``mass_numb`` and exactly one annotated parameter, then ``Z`` may
    be used to provide the |charge number| and ``mass_numb`` may be used
    to provide the mass number. Making ``Z`` and ``mass_numb``
    |keyword-only| reduces the chances of confusion or mistakes.

    >>> @particle_input
    ... def make_particle(particle: ParticleLike, *, Z=None, mass_numb=None):
    ...     return particle
    >>> make_particle("H", Z=0, mass_numb=3)
    Particle("T 0+")

    Instance methods can also be decorated with |particle_input| as long
    as the first argument (representing the instance itself) is named
    ``self`` following standard convention.

    >>> class SampleClass:
    ...
    ...     @particle_input
    ...     def __init__(self, particle: ParticleLike):
    ...         self.particle = particle
    ...
    ...     @particle_input
    ...     def return_particle(self, new_particle: ParticleLike):
    ...         return new_particle

    >>> instance = SampleClass("α")
    >>> instance.particle
    Particle("He-4 2+")
    >>> instance.return_particle("T+")
    Particle("T 1+")

    The ``allow_custom_particles`` and ``allow_particle_lists`` keyword
    arguments indicate whether ``callable_`` should accept
    |CustomParticle| and |ParticleList| objects, respectively.

    >>> @particle_input(allow_custom_particles=False, allow_particle_lists=False)
    ... def get_atomic_number(particle: ParticleLike):
    ...     return particle.atomic_number
    >>> get_atomic_number("Fe")
    26

    The ``require``, ``any_of``, and ``exclude`` keyword arguments may
    be used to specify categorization criteria that a particle must be
    consistent with. For more details, please refer to the docstring for
    `~plasmapy.particles.particle_class.Particle.is_category`.

    >>> @particle_input(require="element", any_of={"charged", "uncharged"})
    ... def return_ionic_level(particle: ParticleLike):
    ...     return particle
    >>> return_ionic_level("Fe-56 0+")
    Particle("Fe-56 0+")

    When the |parameter| is named ``element``, ``isotope``, ``ion``, or
    ``ionic_level``, then the corresponding |argument| must be
    consistent with the name. When the |parameter| is named ``ion`` or
    ``ionic_charge``, then the particle(s) may also be neutral atoms as
    long as the |charge number| is explicitly defined.

    >>> @particle_input
    ... def mass_number(isotope: ParticleLike):
    ...     return isotope.mass_number
    >>> mass_number("D")
    2
    """

    # The following pattern comes from the docs for wrapt, and requires
    # that the arguments to the decorator are keyword-only.

    if callable_ is None:
        return functools.partial(
            particle_input,
            require=require,
            any_of=any_of,
            exclude=exclude,
            allow_custom_particles=allow_custom_particles,
            allow_particle_lists=allow_particle_lists,
        )

    particle_validator = _ParticleInput(
        callable_=callable_,
        require=require,
        any_of=any_of,
        exclude=exclude,
        allow_custom_particles=allow_custom_particles,
        allow_particle_lists=allow_particle_lists,
    )

    @wrapt.decorator
    def wrapper(
        callable__: Callable, instance: Any, args: Tuple, kwargs: Dict[str, Any]
    ):
        new_kwargs = particle_validator.process_arguments(args, kwargs, instance)
        return callable__(**new_kwargs)

    return wrapper(callable_)

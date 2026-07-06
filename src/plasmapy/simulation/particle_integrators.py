"""
Particle movement integrators, for particle simulations.

These do not have `astropy.units` support, choosing instead to
limit overhead and increase performance.

They act in-place on position and velocity arrays to reduce
memory allocation.
"""

__all__ = [
    "AbstractIntegrator",
    "BorisIntegrator",
    "RelativisticBorisIntegrator",
    "RelativisticBorisIntegratorRRF",
]

from abc import ABC, abstractmethod

import astropy.constants as const
import numpy as np

_c = const.c


class AbstractIntegrator(ABC):
    """Outlines the necessary methods to define a particle integrator."""

    @property
    @abstractmethod
    def is_relativistic(self):
        r"""
        Property representing whether an integrator incorporates relativistic
        corrections.
        """
        ...

    @staticmethod
    @abstractmethod
    def push(x, v, B, E, q, m, dt):
        r"""
        The method for applying a push to the specified ensemble of particles.

        The passed parameters should in general not be `~astropy.units.Quantity`
        objects, but rather their SI values. Specific user implementations may vary.
        """
        ...


class BorisIntegrator(AbstractIntegrator):
    """The explicit Boris pusher."""

    @property
    def is_relativistic(self) -> bool:
        r"""
        The explicit Boris pusher is not relativistic.
        """
        return False

    @staticmethod
    def push(x, v, B, E, q, m, dt):
        r"""
        Parameters
        ----------
        x : `~numpy.ndarray`
            Particle position at full timestep, in SI (meter) units.

        v : `~numpy.ndarray`
            Particle velocity at half timestep, in SI (meter/second) units.

        B : `~numpy.ndarray`
            Magnetic field at full timestep, in SI (tesla) units.

        E : `float`
            Electric field at full timestep, in SI (V/m) units.

        q : `float`
            Particle charge, in SI (coulomb) units.

        m : `float`
            Particle mass, in SI (kg) units.

        dt : `float`
            Timestep, in SI (second) units.

        Returns
        -------
        x : `~numpy.ndarray`
            Particle position x after Boris push algorithm, in SI (meter) units.

        v : `~numpy.ndarray`
            Particle velocity after Boris push algorithm, in SI (meter/second) units.

        Examples
        --------
        >>> B = np.array([[0.0, 0.0, 5.0]])
        >>> E = np.array([[0.0, 0.0, 0.0]])
        >>> x_t0 = np.array([[0.0, 0.0, 0.0]])
        >>> v_t0 = np.array([[5.0, 0.0, 0.0]])
        >>> x_t1, v_t1 = BorisIntegrator.push(
        ...     x=x_t0, v=v_t0, B=B, E=E, q=1.0, m=1.0, dt=0.01
        ... )
        >>> x_t1
        array([[ 0.04993754, -0.00249844,  0.        ]])
        >>> v_t1
        array([[ 4.9937539 , -0.24984385,  0.        ]])
        >>> x_t1, v_t1 = BorisIntegrator.push(
        ...     x=x_t0, v=v_t0, B=B, E=E, q=1.0, m=1.0, dt=0.01
        ... )
        >>> x_t1
        array([[ 0.04993754, -0.00249844,  0.        ]])
        >>> v_t1
        array([[ 4.9937539 , -0.24984385,  0.        ]])
        >>> # B parallel to v
        >>> B = np.array([[0.0, 0.0, 5.0]])
        >>> v_t0 = np.array([[0.0, 0.0, 5.0]])
        >>> x_t0 = np.array([[0.0, 0.0, 0.0]])
        >>> x_t1, v_t1 = BorisIntegrator.push(
        ...     x=x_t0, v=v_t0, B=B, E=E, q=1.0, m=1.0, dt=0.01
        ... )
        >>> # no rotation of vector v
        >>> v_t1
        array([[0., 0., 5.]])
        >>> x_t1
        array([[0.  , 0.  , 0.05]])
        >>> # B perpendicular to v
        >>> B = np.array([[5.0, 0.0, 0.0]])
        >>> v_t0 = np.array([[0.0, 5.0, 0.0]])
        >>> x_t0 = np.array([[0.0, 0.0, 0.0]])
        >>> x_t1, v_t1 = BorisIntegrator.push(
        ...     x=x_t0, v=v_t0, B=B, E=E, q=1.0, m=1.0, dt=0.01
        ... )
        >>> # rotation of vector v
        >>> v_t1
        array([[ 0.        ,  4.9937539 , -0.24984385]])
        >>> x_t1
        array([[ 0.        ,  0.04993754, -0.00249844]])
        >>> # nonzero E and zero B
        >>> E = np.array([[1.0, 1.0, 1.0]])
        >>> B = np.array([[0.0, 0.0, 0.0]])
        >>> v_t0 = np.array([[0.0, 0.0, 5.0]])
        >>> x_t0 = np.array([[0.0, 0.0, 0.0]])
        >>> x_t1, v_t1 = BorisIntegrator.push(
        ...     x=x_t0, v=v_t0, B=B, E=E, q=1.0, m=1.0, dt=0.01
        ... )
        >>> v_t1
        array([[0.01, 0.01, 5.01]])
        >>> x_t1
        array([[0.0001, 0.0001, 0.0501]])
        >>> x_t1, v_t1 = BorisIntegrator.push(
        ...     x=x_t0, v=v_t0, B=B, E=E, q=1.0, m=1.0, dt=0.01
        ... )
        >>> v_t1
        array([[0.02, 0.02, 5.02]])
        >>> x_t1
        array([[0.0001, 0.0001, 0.0501]])

        Notes
        -----
        The Boris algorithm :cite:p:`boris:1970` is the standard volume preserving
        algorithm for particle movement in plasma physics. See
        :cite:t:`birdsall:2004` for more details, and this `page on the
        Boris method <https://www.particleincell.com/2011/vxb-rotation>`__
        for a nice overview.

        Conceptually, the algorithm has three phases:

        1. Add half the impulse from electric field.
        2. Rotate the particle velocity about the direction of the magnetic
           field.
        3. Add the second half of the impulse from the electric field.

        This ends up causing the magnetic field action to be properly "centered" in
        time, and the algorithm being volume preserving ensures that error in energy
        remains bounded. See this paper <https://pubs.aip.org/aip/pop/article/20/8/084503/317652/Why-is-Boris-algorithm-so-good>
        for more detail.
        """
        hqmdt = 0.5 * dt * q / m
        vminus = v + hqmdt * E

        # rotate to add magnetic field
        t = B * hqmdt
        s = 2 * t / (1 + (t * t).sum(axis=1, keepdims=True))
        vprime = vminus + np.cross(vminus, t)
        vplus = vminus + np.cross(vprime, s)

        # add second half of electric impulse
        v = vplus + hqmdt * E
        return x + v * dt, v


class RelativisticBorisIntegrator(AbstractIntegrator):
    """The explicit Boris pusher, including relativistic corrections."""

    @property
    def is_relativistic(self) -> bool:
        r"""
        The push implementation of the Boris push algorithm includes relativistic
        corrections.
        """
        return True

    @staticmethod
    def _boris_lorentz_push(v, B, E, q, m, dt):
        r"""
        Helper method applying the relativistic Boris push to the proper velocity.

        This method carries out the Lorentz-force part of the Boris push. It does not convert
        back to a coordinate velocity or advance position, each push() method completes this
        function. Both proper velocities are returned so the RRF pusher can estimate total
        momentum & velocity p{n} and v{n} and apply the additional momentum corrections. The
        regular relativistic boris pusher can be called to complete the timestep from uvel_L
        alone, and discards the pre Lorentz proper velocity (uvel)

        Parameters
        ----------
        v : `~numpy.ndarray`
            particle velocity at half timestep, in SI (meter/second) units.
        B : `~numpy.ndarray`
            magnetic field at full timestep, in SI (tesla) units.
        E : float
            electric field at full timestep, in SI (V/m) units.
        q : float
            particle charge, in SI (Coulomb) units.
        m : float
            particle mass, in SI (kg) units.
        dt : float
            timestep, in SI (second) units.

        Returns
        -------
        uvel : `~numpy.ndarray`
            proper velocity before the push, p{n-1/2}/m, in SI (meter/second)
            units.
        uvel_L : `~numpy.ndarray`
            proper velocity after the Lorentz-force push, p_L{n+1/2}/m, in SI
            (meter/second) units.
        """
        # need float value for c
        c = _c.si.value
        γ = 1 / np.sqrt(1 - (np.linalg.norm(v, axis=1, keepdims=True) / c) ** 2)
        uvel = v * γ
        uvel_minus = uvel + q * E * dt / (2 * m)
        γ1 = np.sqrt(1 + (np.linalg.norm(uvel_minus, axis=1, keepdims=True) / c) ** 2)
        t = q * B * dt / (2 * γ1 * m)
        s = 2 * t / (1 + (t * t).sum(axis=1, keepdims=True))
        uvel_prime = uvel_minus + np.cross(uvel_minus, t)
        uvel_plus = uvel_minus + np.cross(uvel_prime, s)
        uvel_L = uvel_plus + q * E * dt / (2 * m)
        # uvel = p{n-1/2}/m is returned with uvel_L = p_L{n+1/2}/m so the RRF
        # pusher can estimate integer-step quantities
        return uvel, uvel_L

    @staticmethod
    def push(x, v, B, E, q, m, dt):
        r"""
        Parameters
        ----------
        x : `~numpy.ndarray`
            particle position at full timestep, in SI (meter) units.
        v : `~numpy.ndarray`
            particle velocity at half timestep, in SI (meter/second) units.
        B : `~numpy.ndarray`
            magnetic field at full timestep, in SI (tesla) units.
        E : float
            electric field at full timestep, in SI (V/m) units.
        q : float
            particle charge, in SI (Coulomb) units.
        m : float
            particle mass, in SI (kg) units.
        dt : float
            timestep, in SI (second) units.

        Notes
        -----
        For the basic overview of this algorithm, see `BorisIntegrator`. This
        version, based on [1]_, applies relativistic corrections such as
        the proper velocity and proper time transformations. These corrections
        affect the leapfrog scheme timestep. The addition of the impulse from
        the electric field must now account for the relativistic mass, resulting
        in an impulse reduced by a factor of :math:`\gamma` relative to the classical
        result. A similar adjustment occurs in the rotation due to the magnetic
        field, as the angle of rotation is roughly inversely proportional to the
        mass of the charged particle. More precisely, it is proportional to :math:`\arctan{(1 / m)}`.
        This means the angle of rotation will also decrease approximately by a
        factor of :math:`\gamma`.

        Keep in mind that the non-relativistic version will be slightly
        faster if you don't encounter velocities in relativistic regimes.

        References
        ----------
        .. [1] C. K. Birdsall, A. B. Langdon, "Plasma Physics via Computer
               Simulation", 2004, p. 58-63
        """
        # discard uvel with _ pre-push proper velocity is only needed by the RRF pusher
        _, uvel_L = RelativisticBorisIntegrator._boris_lorentz_push(v, B, E, q, m, dt)
        # You can show that this expression is equivalent to calculating
        # v_new  then calculating γnew using the usual formula
        γ2 = np.sqrt(
            1 + (np.linalg.norm(uvel_L, axis=1, keepdims=True) / _c.si.value) ** 2,
        )

        v = uvel_L / γ2

        return x + v * dt, v


class RelativisticBorisIntegratorRRF(RelativisticBorisIntegrator):
    """Relativistic Boris Pusher with Tamburini Radiation Reaction Force Implementation."""

    @staticmethod
    def push(x, v, B, E, q, m, dt):
        r"""
        Parameters
        ----------
        x : `~numpy.ndarray`
            particle position at full timestep, in SI (meter) units.
        v : `~numpy.ndarray`
            particle velocity at half timestep, in SI (meter/second) units.
        B : `~numpy.ndarray`
            magnetic field at full timestep, in SI (tesla) units.
        E : float
            electric field at full timestep, in SI (V/m) units.
        q : float
            particle charge, in SI (Coulomb) units.
        m : float
            particle mass, in SI (kg) units.
        dt : float
            timestep, in SI (second) units.

        Returns
        -------
        x : `~numpy.ndarray`
            Particle position x after Boris push algorithm, in SI (meter) units.
        v : `~numpy.ndarray`
            Particle velocity after Boris push algorithm, in SI (meter/second) units.

        Notes
        -----
        This method is the specific push including a RRF damping effect within the
        RelativisticBorisIntegratorRRF Class. This method takes into account a
        radiation reaction force under various important assumptions. The process of
        including this radiative damping force leaves the standard Boris algorithm
        unchanged and works in the following steps.

        1. The ``_boris_lorentz_push`` helper is executed; the Lorentz-pushed
           momentum p_L{n+1/2} and the initial momentum p{n-1/2} are evaluated
           and stored.
        2. The radiation reaction force is then estimated at the n integer
           step state utilizing these values and various assumptions.
        3. This RRF is then applied as a momentum kick with the following
           equation serving as the basis: p{n+1/2} = p_L{n+1/2} + f_R{n} * dt.

        This RRF implementation follows the post-Boris scheme of
        :cite:t:`tamburini:2010`.

        This paper includes many important stipulations and utilizes the
        Landau-Lifshitz radiation reaction approximation model. :cite:t:`landau:1975`
        The LL model simplifies the Abraham-Lorentz-Dirac (ALD) equation, which contains a
        third-order time derivative (da/dt, or a jerk term) that produces clearly
        unphysical solutions such as runaway (infinitely accelerating), and
        pre-acceleration solutions. This LL model eliminates the jerk term by
        substituting the instantaneous acceleration due to the Lorentz Force,
        which can only be approximated under certain conditions.

        1. "The acceleration of particles is dominated by the Lorentz force,
           with the RR force giving a smaller, albeit non negligible
           contribution." F_L >> F_R :cite:p:`tamburini:2010`

        2. The rest-frame field must stay well below the Schwinger critical
           field (E_cr ~ 1.3e18 V/m, B_cr ~ 4.4e9 T)
           :cite:p:`tamburini:2010` :cite:p:`landau:1975`

        3. Neglected f_R terms (Tamburini Eq. 6): the field-derivative term
           is dropped and the electron spin force is likewise neglected.
           :cite:p:`tamburini:2010`

        Examples
        --------
        >>> x = np.array([[0.0, 0.0, 0.0]])
        >>> v = np.array([[1.0e8, 0.0, 0.0]])
        >>> B = np.array([[0.0, 0.0, 1.0e3]])
        >>> E = np.zeros((1, 3))
        >>> x_new, v_new = RelativisticBorisIntegratorRRF.push(
        ...     x, v, B, E, q=-1.6e-19, m=9.1e-31, dt=1.0e-13
        ... )
        >>> bool(np.linalg.norm(v_new) <= np.linalg.norm(v))  # radiation can't add velo
        True
        """
        # need float value for c
        c = _c.si.value
        # relativistic boris push method to determine Lorentz momentum boost p_L{n+1/2}
        # and total momentum @ n - 1/2:    p{n-1/2}
        uvel, uvel_L = RelativisticBorisIntegratorRRF._boris_lorentz_push(
            v, B, E, q, m, dt
        )
        # proper velocity obtained w/ out RR
        # at this point we have same result as RelativisticBorisIntegrator.push(),
        # just storing value and noting: uvel = p{n-1/2}/m  ;  uvel_L = p_L{n+1/2}/m
        # use these stored p{n-1/2}/m and p_L{n+1/2}/m values to estimate total
        # momentum & velocity p{n} and v{n} at the integer step n (from eqn --> [12])
        # p{n} =~ .5(p_L{n+1/2} + p{n-1/2})  [12]
        u_n = 1 / 2 * (uvel_L + uvel)
        # gamma calculation
        γn = np.sqrt(
            1 + (np.linalg.norm(u_n, axis=1, keepdims=True) / c) ** 2,
        )
        # integer step approximate velocity: v{n} =~ p{n}/γ{n} [12]
        v_n = u_n / γn
        # take your p{n} and subsequently v{n} value to calculate RRF @ integer step n
        f_R = RelativisticBorisIntegratorRRF.rrf_full(v_n, B, E, q, m)
        # finally utilizing this new calculated RRF, can apply post-boris kick like described in eqn [7-11]
        # (p{n+1/2} - p{n-1/2}) / dt  = f{n} = f_L{n} + f_R{n} --> p{n+1/2} = p_L{n+1/2} + f_R{n} * dt
        uvel_final = uvel_L + f_R * dt / m
        # gamma calculation
        γ2 = np.sqrt(
            1 + (np.linalg.norm(uvel_final, axis=1, keepdims=True) / c) ** 2,
        )
        # convert to velo and advance position from this value as described in final eqn above
        v_final = uvel_final / γ2
        return x + v_final * dt, v_final

    @staticmethod
    def rrf_full(v, B, E, q, m):
        r"""
        Radiation-reaction force from the Landau-Lifshitz approximation,
        following :cite:t:`tamburini:2010`.

        Parameters
        ----------
        v : `~numpy.ndarray`
            particle velocity at the integer timestep, in SI (meter/second)
            units.
        B : `~numpy.ndarray`
            magnetic field at full timestep, in SI (tesla) units.
        E : float
            electric field at full timestep, in SI (V/m) units.
        q : float
            particle charge, in SI (Coulomb) units.
        m : float
            particle mass, in SI (kg) units.

        Returns
        -------
        f_R : `~numpy.ndarray`
            radiation-reaction force at the integer timestep, in SI (newton)
            units.

        Examples
        --------
        >>> v = np.array([[1.0e7, 0.0, 0.0]])
        >>> B = np.array([[0.0, 0.0, 1.0]])
        >>> E = np.zeros((1, 3))
        >>> f_R = RelativisticBorisIntegratorRRF.rrf_full(
        ...     v, B, E, q=-1.6e-19, m=9.1e-31
        ... )
        >>> f_R.shape
        (1, 3)
        >>> bool(np.dot(f_R[0], v[0]) < 0)  # the force opposes the velocity (a drag)
        True
        """
        # need float value for c
        c = _c.si.value
        # gamma calculation [[1]]
        γ = 1 / np.sqrt(
            1 - (np.linalg.norm(v, axis=1, keepdims=True) / c) ** 2,
        )
        # constant out front [[C^2 s / kg]]
        k = q**4 / (6 * np.pi * const.eps0.si.value * m**2 * c**3)
        # Lorentz Force  f_L ≡ -(E+v×B) [[V/m]] ==> [[N/C]]
        f_L = -(E + np.cross(v, B))
        # Lorentz Force squared (helpful quantity from Tamburini et al.) [[V^2/m^2]]
        f_L_squared = (f_L * f_L).sum(axis=1, keepdims=True)
        # (v dot E) (another helpful quantity from Tamburini et al.) [[V/s]]
        v_dotproduct_E = (v * E).sum(axis=1, keepdims=True)
        # term 1 eqn (6) [[kg^2 m / (s^3 C^2)]]
        rrf_term1 = np.cross(f_L, B) - (v_dotproduct_E / c**2) * E
        # term 2 eqn (6) [[kg^2 m / (s^3 C^2)]]
        rrf_term2 = (γ**2 / c**2) * (f_L_squared - v_dotproduct_E**2 / c**2) * v
        # final f_R [[C^2 s/kg]][[kg^2 m / (s^3 C^2)]] --> [[N]]
        return -k * (rrf_term1 + rrf_term2)

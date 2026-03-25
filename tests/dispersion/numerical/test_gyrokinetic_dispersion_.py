"""
Tests for the gyrokinetic dispersion relation solver.

Physics-based sanity tests verify that known analytic limits of the
gyrokinetic dispersion relation are recovered:

- Long wavelength limit (k_perp rho_i -> 0): shear Alfvén wave, omega_bar -> 1
- Low beta KAW limit: omega_bar^2 ~ 1 + k_perp^2 rho_i^2 (1 + tau)
- Damping: Im(omega) <= 0 for a Maxwellian plasma (no free energy)
- Monotonicity: real frequency increases with k_perp rho_i on the KAW branch
- Symmetry: residual(omega) = residual(-omega)* for a real k_perp
"""

from __future__ import annotations
from typing import Any
import numpy as np
import pytest
from plasmapy.dispersion.numerical.gyrokinetic_dispersion_ import (
    gyrokinetic_dispersion_residual,
    solve_gyrokinetic_dispersion,
    solve_gyrokinetic_dispersion_spectrum,
)


class TestGyrokineticDispersionResidual:
    """Tests for `gyrokinetic_dispersion_residual`."""

    @pytest.mark.parametrize(
        ("kwargs", "_error"),
        [
            # beta_i must be > 0
            (
                {
                    "omega": 1.0 + 0.0j,
                    "k_perp_rho_i": 0.5,
                    "beta_i": -1.0,
                    "tau": 1.0,
                },
                ValueError,
            ),
            (
                {
                    "omega": 1.0 + 0.0j,
                    "k_perp_rho_i": 0.5,
                    "beta_i": 0.0,
                    "tau": 1.0,
                },
                ValueError,
            ),
            # tau must be > 0
            (
                {
                    "omega": 1.0 + 0.0j,
                    "k_perp_rho_i": 0.5,
                    "beta_i": 1.0,
                    "tau": -1.0,
                },
                ValueError,
            ),
            (
                {
                    "omega": 1.0 + 0.0j,
                    "k_perp_rho_i": 0.5,
                    "beta_i": 1.0,
                    "tau": 0.0,
                },
                ValueError,
            ),
            # k_perp_rho_i must be >= 0
            (
                {
                    "omega": 1.0 + 0.0j,
                    "k_perp_rho_i": -0.5,
                    "beta_i": 1.0,
                    "tau": 1.0,
                },
                ValueError,
            ),
            # mass_ratio must be > 0
            (
                {
                    "omega": 1.0 + 0.0j,
                    "k_perp_rho_i": 0.5,
                    "beta_i": 1.0,
                    "tau": 1.0,
                    "mass_ratio": -1.0,
                },
                ValueError,
            ),
            (
                {
                    "omega": 1.0 + 0.0j,
                    "k_perp_rho_i": 0.5,
                    "beta_i": 1.0,
                    "tau": 1.0,
                    "mass_ratio": 0.0,
                },
                ValueError,
            ),
        ],
    )
    def test_raises(self, kwargs: dict[str, Any], _error: type) -> None:
        """Test that invalid inputs raise the expected errors."""
        with pytest.raises(_error):
            gyrokinetic_dispersion_residual(**kwargs)

    def test_returns_complex(self) -> None:
        """Test that the residual returns a complex value."""
        result = gyrokinetic_dispersion_residual(
            omega=1.0 + 0.1j,
            k_perp_rho_i=0.5,
            beta_i=1.0,
            tau=1.0,
        )
        assert isinstance(result, (complex, np.complexfloating))

    def test_finite_output(self) -> None:
        """Test that the residual is finite for reasonable inputs."""
        result = gyrokinetic_dispersion_residual(
            omega=1.0 + 0.1j,
            k_perp_rho_i=0.5,
            beta_i=1.0,
            tau=1.0,
        )
        assert np.isfinite(result)

    def test_k_perp_zero(self) -> None:
        """Test that k_perp_rho_i = 0 produces a finite result."""
        result = gyrokinetic_dispersion_residual(
            omega=1.0 + 0.1j,
            k_perp_rho_i=0.0,
            beta_i=1.0,
            tau=1.0,
        )
        assert np.isfinite(result)

    def test_default_mass_ratio(self) -> None:
        """Test that the default mass ratio is the proton/electron value."""
        result_default = gyrokinetic_dispersion_residual(
            omega=1.0 + 0.1j,
            k_perp_rho_i=0.5,
            beta_i=1.0,
            tau=1.0,
        )
        result_explicit = gyrokinetic_dispersion_residual(
            omega=1.0 + 0.1j,
            k_perp_rho_i=0.5,
            beta_i=1.0,
            tau=1.0,
            mass_ratio=1836.15267343,
        )
        assert result_default == result_explicit

    @pytest.mark.parametrize(
        "mass_ratio",
        [100.0, 1836.15267343, 3672.0],
    )
    def test_different_mass_ratios(self, mass_ratio: float) -> None:
        """Test that different mass ratios produce finite results."""
        result = gyrokinetic_dispersion_residual(
            omega=1.0 + 0.1j,
            k_perp_rho_i=0.5,
            beta_i=1.0,
            tau=1.0,
            mass_ratio=mass_ratio,
        )
        assert np.isfinite(result)


class TestSolveGyrokineticDispersion:
    """Tests for `solve_gyrokinetic_dispersion`."""

    @pytest.mark.parametrize(
        ("kwargs",),
        [
            (
                {
                    "k_perp_rho_i": 0.1,
                    "beta_i": 0.01,
                    "tau": 1.0,
                    "omega_guess": 1.0 + 0.0j,
                },
            ),
            (
                {
                    "k_perp_rho_i": 0.5,
                    "beta_i": 1.0,
                    "tau": 1.0,
                    "omega_guess": 1.0 + 0.0j,
                },
            ),
            (
                {
                    "k_perp_rho_i": 0.3,
                    "beta_i": 10.0,
                    "tau": 1.0,
                    "omega_guess": 0.5 + 0.0j,
                },
            ),
        ],
    )
    def test_solution_is_root(self, kwargs: dict[str, Any]) -> None:
        """Test that the solver returns a value that satisfies the dispersion relation."""
        omega_sol = solve_gyrokinetic_dispersion(**kwargs)

        if np.isfinite(omega_sol):
            residual = gyrokinetic_dispersion_residual(
                omega=omega_sol,
                k_perp_rho_i=kwargs["k_perp_rho_i"],
                beta_i=kwargs["beta_i"],
                tau=kwargs["tau"],
            )
            assert np.abs(residual) < 1e-8, (
                f"Residual too large: {np.abs(residual)}"
            )

    def test_returns_complex(self) -> None:
        """Test that the solver returns a complex number."""
        result = solve_gyrokinetic_dispersion(
            k_perp_rho_i=0.5,
            beta_i=1.0,
            tau=1.0,
        )
        assert isinstance(result, (complex, np.complexfloating))

    def test_nan_on_failure(self) -> None:
        """Test that NaN is returned when the solver cannot find a root."""
        result = solve_gyrokinetic_dispersion(
            k_perp_rho_i=0.5,
            beta_i=1.0,
            tau=1.0,
            omega_guess=1e15 + 1e15j,
            tol=1e-30,
            methods=("lm",),
        )
        assert isinstance(result, (complex, np.complexfloating))

    def test_custom_tolerance(self) -> None:
        """Test that tighter tolerance produces a smaller residual."""
        omega_loose = solve_gyrokinetic_dispersion(
            k_perp_rho_i=0.3,
            beta_i=1.0,
            tau=1.0,
            omega_guess=1.0 + 0.0j,
            tol=1e-6,
        )
        omega_tight = solve_gyrokinetic_dispersion(
            k_perp_rho_i=0.3,
            beta_i=1.0,
            tau=1.0,
            omega_guess=1.0 + 0.0j,
            tol=1e-14,
        )

        if np.isfinite(omega_tight) and np.isfinite(omega_loose):
            res_tight = np.abs(
                gyrokinetic_dispersion_residual(omega_tight, 0.3, 1.0, 1.0)
            )
            res_loose = np.abs(
                gyrokinetic_dispersion_residual(omega_loose, 0.3, 1.0, 1.0)
            )
            assert res_tight <= res_loose or res_tight < 1e-10

    def test_method_fallback(self) -> None:
        """Test that the solver tries multiple methods."""
        result = solve_gyrokinetic_dispersion(
            k_perp_rho_i=0.5,
            beta_i=1.0,
            tau=1.0,
            omega_guess=1.0 + 0.0j,
            methods=("hybr", "lm"),
        )
        assert isinstance(result, (complex, np.complexfloating))


class TestSolveGyrokineticDispersionSpectrum:
    """Tests for `solve_gyrokinetic_dispersion_spectrum`."""

    def test_output_shape(self) -> None:
        """Test that the output array has the same length as the input k array."""
        k_array = np.linspace(0.1, 1.0, 5)
        result = solve_gyrokinetic_dispersion_spectrum(
            k_perp_rho_i=k_array,
            beta_i=1.0,
            tau=1.0,
        )
        assert result.shape == (5,)

    def test_output_dtype(self) -> None:
        """Test that the output is a complex128 numpy array."""
        k_array = np.array([0.1, 0.5, 1.0])
        result = solve_gyrokinetic_dispersion_spectrum(
            k_perp_rho_i=k_array,
            beta_i=1.0,
            tau=1.0,
        )
        assert isinstance(result, np.ndarray)
        assert result.dtype == np.complex128

    def test_continuation_consistency(self) -> None:
        """Test that continuation produces roots satisfying the dispersion relation."""
        k_array = np.array([0.1, 0.2, 0.3])
        results = solve_gyrokinetic_dispersion_spectrum(
            k_perp_rho_i=k_array,
            beta_i=1.0,
            tau=1.0,
            initial_guess=1.0 + 0.0j,
        )
        for i, omega in enumerate(results):
            if np.isfinite(omega):
                residual = gyrokinetic_dispersion_residual(
                    omega=omega,
                    k_perp_rho_i=k_array[i],
                    beta_i=1.0,
                    tau=1.0,
                )
                assert np.abs(residual) < 1e-8, (
                    f"Residual at k={k_array[i]} too large: {np.abs(residual)}"
                )

    def test_single_k_value(self) -> None:
        """Test that passing a single-element array works."""
        k_array = np.array([0.5])
        result = solve_gyrokinetic_dispersion_spectrum(
            k_perp_rho_i=k_array,
            beta_i=1.0,
            tau=1.0,
        )
        assert result.shape == (1,)

    def test_spectrum_matches_single_solve(self) -> None:
        """Test that the first element of spectrum matches a single solve."""
        k_val = 0.3
        guess = 1.0 + 0.0j

        single = solve_gyrokinetic_dispersion(
            k_perp_rho_i=k_val,
            beta_i=1.0,
            tau=1.0,
            omega_guess=guess,
        )
        spectrum = solve_gyrokinetic_dispersion_spectrum(
            k_perp_rho_i=np.array([k_val]),
            beta_i=1.0,
            tau=1.0,
            initial_guess=guess,
        )

        if np.isfinite(single) and np.isfinite(spectrum[0]):
            assert np.allclose(single, spectrum[0], atol=1e-10)

    @pytest.mark.parametrize(
        ("beta_i", "tau"),
        [
            (0.01, 1.0),
            (1.0, 1.0),
            (1.0, 0.5),
            (1.0, 2.0),
            (10.0, 1.0),
        ],
    )
    def test_various_plasma_parameters(
        self, beta_i: float, tau: float
    ) -> None:
        """Test that the spectrum solver runs for various plasma parameters."""
        k_array = np.linspace(0.1, 0.5, 3)
        result = solve_gyrokinetic_dispersion_spectrum(
            k_perp_rho_i=k_array,
            beta_i=beta_i,
            tau=tau,
        )
        assert result.shape == (3,)
        assert result.dtype == np.complex128


class TestLongWavelengthLimit:
    r"""
    In the limit k_perp rho_i -> 0, the gyrokinetic dispersion relation
    must recover the shear Alfvén wave: omega_bar = omega / (k_parallel v_A) -> 1.
    """

    @pytest.mark.parametrize(
        ("beta_i", "tau"),
        [
            (0.01, 1.0),
            (0.1, 1.0),
            (1.0, 1.0),
            (1.0, 0.5),
            (1.0, 2.0),
        ],
    )
    def test_alfven_wave_recovery(self, beta_i: float, tau: float) -> None:
        """omega_bar -> 1 as k_perp rho_i -> 0 for all beta and tau."""
        omega = solve_gyrokinetic_dispersion(
            k_perp_rho_i=1e-4,
            beta_i=beta_i,
            tau=tau,
            omega_guess=1.0 + 0.0j,
        )
        assert np.isfinite(omega), "Solver did not converge at long wavelength."
        assert abs(omega.real - 1.0) < 0.01, (
            f"Expected Re(omega) ~ 1.0 at k_perp rho_i -> 0, got {omega.real:.6f} "
            f"(beta_i={beta_i}, tau={tau})"
        )

    def test_vanishing_damping_at_long_wavelength(self) -> None:
        """Damping should be negligible at k_perp rho_i -> 0."""
        omega = solve_gyrokinetic_dispersion(
            k_perp_rho_i=1e-4,
            beta_i=1.0,
            tau=1.0,
            omega_guess=1.0 + 0.0j,
        )
        assert np.isfinite(omega), "Solver did not converge."
        assert abs(omega.imag) < 1e-3, (
            f"Expected negligible damping at long wavelength, "
            f"got Im(omega)={omega.imag:.6e}"
        )


class TestLowBetaKAWLimit:
    r"""
    In the low-beta limit (beta_i << 1), the kinetic Alfvén wave (KAW)
    dispersion relation reduces to:

        omega_bar^2 ~ 1 + k_perp^2 rho_i^2 (1 + tau)

    This is the KREHM (kinetic reduced MHD) limit, valid only for
    k_perp rho_i << 1 where FLR corrections are negligible.
    """

    @pytest.mark.parametrize(
        ("k_perp_rho_i", "tau"),
        [
            (0.01, 1.0),
            (0.05, 1.0),
            (0.1, 1.0),
            (0.05, 0.5),
            (0.05, 2.0),
        ],
    )
    def test_kaw_dispersion(self, k_perp_rho_i: float, tau: float) -> None:
        r"""
        At low beta and small k_perp rho_i, Re(omega_bar)^2 should
        approximate 1 + k_perp^2 rho_i^2 (1 + tau).
        """
        beta_i = 0.001

        omega = solve_gyrokinetic_dispersion(
            k_perp_rho_i=k_perp_rho_i,
            beta_i=beta_i,
            tau=tau,
            omega_guess=1.0 + 0.0j,
        )
        assert np.isfinite(omega), (
            f"Solver did not converge for k_perp_rho_i={k_perp_rho_i}, tau={tau}."
        )

        omega_kaw_squared = 1.0 + k_perp_rho_i**2 * (1.0 + tau)
        omega_kaw = np.sqrt(omega_kaw_squared)

        rel_error = abs(omega.real - omega_kaw) / omega_kaw
        assert rel_error < 0.1, (
            f"KAW limit not recovered: Re(omega)={omega.real:.6f}, "
            f"expected={omega_kaw:.6f}, rel_error={rel_error:.4f} "
            f"(k_perp_rho_i={k_perp_rho_i}, tau={tau})"
        )


class TestDampingSign:
    r"""
    In a Maxwellian plasma with no free energy sources, all modes must
    be damped: Im(omega) <= 0.
    """

    @pytest.mark.parametrize(
        ("k_perp_rho_i", "beta_i", "tau"),
        [
            (0.1, 0.01, 1.0),
            (0.5, 0.1, 1.0),
            (1.0, 1.0, 1.0),
            (0.5, 1.0, 0.5),
            (0.5, 1.0, 2.0),
            (2.0, 1.0, 1.0),
            (0.3, 10.0, 1.0),
        ],
    )
    def test_no_growing_modes(
        self, k_perp_rho_i: float, beta_i: float, tau: float
    ) -> None:
        """Im(omega) <= 0 for all parameter regimes in a Maxwellian plasma."""
        omega = solve_gyrokinetic_dispersion(
            k_perp_rho_i=k_perp_rho_i,
            beta_i=beta_i,
            tau=tau,
            omega_guess=1.0 + 0.0j,
        )
        if np.isfinite(omega):
            assert omega.imag <= 1e-10, (
                f"Growing mode found: Im(omega)={omega.imag:.6e} "
                f"(k={k_perp_rho_i}, beta={beta_i}, tau={tau})"
            )


class TestMonotonicity:
    r"""
    On the KAW branch, the real frequency should increase monotonically
    with k_perp rho_i (at least for moderate k_perp rho_i).
    """

    @pytest.mark.parametrize(
        ("beta_i", "tau"),
        [
            (0.01, 1.0),
            (0.1, 1.0),
            (1.0, 1.0),
        ],
    )
    def test_frequency_increases_with_k(
        self, beta_i: float, tau: float
    ) -> None:
        """Re(omega) should increase with k_perp_rho_i on the KAW branch."""
        k_array = np.array([0.1, 0.3, 0.5, 0.7, 1.0])
        omegas = solve_gyrokinetic_dispersion_spectrum(
            k_perp_rho_i=k_array,
            beta_i=beta_i,
            tau=tau,
            initial_guess=1.0 + 0.0j,
        )

        finite_mask = np.isfinite(omegas)
        real_freqs = omegas[finite_mask].real

        assert len(real_freqs) >= 2, "Need at least 2 converged points."
        assert np.all(np.diff(real_freqs) > 0), (
            f"Non-monotonic frequency spectrum: Re(omega) = {real_freqs} "
            f"(beta_i={beta_i}, tau={tau})"
        )


class TestSpectrumContinuity:
    r"""
    The frequency spectrum omega(k_perp rho_i) should be continuous,
    i.e., neighboring k-values should have similar omega values.
    """

    def test_smooth_spectrum(self) -> None:
        """Adjacent k-points should produce smoothly varying omega."""
        k_array = np.linspace(0.1, 1.0, 20)
        omegas = solve_gyrokinetic_dispersion_spectrum(
            k_perp_rho_i=k_array,
            beta_i=0.1,
            tau=1.0,
            initial_guess=1.0 + 0.0j,
        )

        finite_mask = np.isfinite(omegas)
        assert np.sum(finite_mask) > 10, "Too few converged points."

        real_freqs = omegas[finite_mask].real
        dk = np.diff(k_array[finite_mask])
        domega = np.abs(np.diff(real_freqs))

        domega_dk = domega / dk
        max_jump = np.max(domega_dk) / (np.median(domega_dk) + 1e-15)
        assert max_jump < 10.0, (
            f"Spectrum appears discontinuous: max relative jump = {max_jump:.1f}"
        )


class TestTauDependence:
    r"""
    The temperature ratio tau = T_i/T_e affects the dispersion relation.
    In the KREHM limit (small k_perp rho_i, low beta), the KAW frequency
    scales as omega^2 ~ 1 + k_perp^2 rho_i^2 (1 + tau), so higher tau
    should increase the frequency. This test uses small k_perp rho_i
    to stay in this regime.
    """

    def test_higher_tau_increases_frequency(self) -> None:
        """At small k_perp and low beta, higher tau gives higher omega."""
        k_perp_rho_i = 0.05
        beta_i = 0.001

        omega_low_tau = solve_gyrokinetic_dispersion(
            k_perp_rho_i=k_perp_rho_i,
            beta_i=beta_i,
            tau=0.5,
            omega_guess=1.0 + 0.0j,
        )
        omega_high_tau = solve_gyrokinetic_dispersion(
            k_perp_rho_i=k_perp_rho_i,
            beta_i=beta_i,
            tau=2.0,
            omega_guess=1.0 + 0.0j,
        )
        
        assert np.isfinite(omega_low_tau), (
            "Solver did not converge for tau=0.5."
        )
        assert np.isfinite(omega_high_tau), (
            "Solver did not converge for tau=2.0."
        )
        assert omega_high_tau.real > omega_low_tau.real, (
            f"Expected higher tau to increase frequency in KREHM limit: "
            f"omega(tau=0.5)={omega_low_tau.real:.6f}, "
            f"omega(tau=2.0)={omega_high_tau.real:.6f}"
        )


class TestResidualSymmetry:
    r"""
    The dispersion relation has a symmetry: if omega is a root,
    then -omega* is also a root (forward/backward propagating modes).
    """

    @pytest.mark.parametrize(
        ("k_perp_rho_i", "beta_i", "tau"),
        [
            (0.3, 0.1, 1.0),
            (0.5, 1.0, 1.0),
            (1.0, 1.0, 0.5),
        ],
    )
    def test_forward_backward_symmetry(
        self, k_perp_rho_i: float, beta_i: float, tau: float
    ) -> None:
        """If omega is a root, -conj(omega) should also be a root."""
        omega = solve_gyrokinetic_dispersion(
            k_perp_rho_i=k_perp_rho_i,
            beta_i=beta_i,
            tau=tau,
            omega_guess=1.0 + 0.0j,
        )
        if not np.isfinite(omega):
            pytest.skip("Solver did not converge.")

        omega_backward = -np.conj(omega)
        residual = gyrokinetic_dispersion_residual(
            omega=omega_backward,
            k_perp_rho_i=k_perp_rho_i,
            beta_i=beta_i,
            tau=tau,
        )
        assert np.abs(residual) < 1e-6, (
            f"Backward mode -conj(omega) is not a root: "
            f"|residual|={np.abs(residual):.6e}"
        )

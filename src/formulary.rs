use numpy::{PyArrayDyn, PyReadonlyArrayDyn, ToPyArray};
use pyo3::prelude::*;
use pyo3::types::PyBool;

#[pyfunction]
pub fn plasma_frequency_lite<'py>(
    py: pyo3::Python<'py>,
    n: PyReadonlyArrayDyn<'_, f64>,
    mass: PyReadonlyArrayDyn<'_, f64>,
    z_mean: PyReadonlyArrayDyn<'_, f64>,
    to_hz: &PyBool,
) -> &'py PyArrayDyn<f64> {
    let n = n.as_array();
    let mass = mass.as_array();
    let z_mean = z_mean.as_array();

    let radicand = &n / (physical_constants::VACUUM_ELECTRIC_PERMITTIVITY * &mass);
    let square_root = radicand.map(|v| v.sqrt());

    let omega_p = &z_mean * physical_constants::ELEMENTARY_CHARGE * &square_root;

    let result = if !to_hz.is_true() {
        omega_p
    } else {
        (&omega_p) / (2.0 * std::f64::consts::PI)
    };

    result.to_pyarray(py)
}

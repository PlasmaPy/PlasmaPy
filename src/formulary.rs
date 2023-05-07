use ndarray::{Array, ArrayView, IxDyn};
use numpy::ndarray::Zip;
use numpy::{PyArrayDyn, PyReadonlyArrayDyn, ToPyArray};
use pyo3::prelude::*;

#[pyfunction]
pub fn plasma_frequency_lite_wrapper<'py>(
    py: Python<'py>,
    n: PyReadonlyArrayDyn<'_, f64>,
    mass: PyReadonlyArrayDyn<'_, f64>,
    z_mean: PyReadonlyArrayDyn<'_, f64>,
    to_hz: PyReadonlyArrayDyn<'_, bool>,
) -> &'py PyArrayDyn<f64> {
    let n = n.as_array();
    let mass = mass.as_array();
    let z_mean = z_mean.as_array();
    let to_hz = to_hz.as_array();

    let omega_p = py.allow_threads(|| plasma_frequency_lite(n, mass, z_mean, to_hz));

    omega_p.to_pyarray(py)
}

fn plasma_frequency_lite(
    n: ArrayView<f64, IxDyn>,
    mass: ArrayView<f64, IxDyn>,
    z_mean: ArrayView<f64, IxDyn>,
    _to_hz: ArrayView<bool, IxDyn>,
) -> Array<f64, IxDyn> {
    let mut n = n.to_owned();

    Zip::from(&mut n)
        .and(&mass)
        .and(&z_mean)
        .par_for_each(|x, &y, &z| {
            *x = z
                * physical_constants::ELEMENTARY_CHARGE
                * (*x / (physical_constants::VACUUM_ELECTRIC_PERMITTIVITY * y)).sqrt()
        });

    n
}

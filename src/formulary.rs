use pyo3::{PyResult, pyfunction};

#[pyfunction]
pub fn my_function() -> PyResult<f32>{
    Ok(21.0)
}

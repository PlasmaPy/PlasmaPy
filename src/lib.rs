use pyo3::prelude::*;

mod formulary;

#[pymodule]
fn rust(_py: Python, module: &PyModule) -> PyResult<()> {
    register_child_module(_py, module)?;
    Ok(())
}

fn register_child_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()> {
    let formulary_module = PyModule::new(py, "formulary")?;
    formulary_module.add_function(wrap_pyfunction!(formulary::my_function, formulary_module)?)?;
    parent_module.add_submodule(formulary_module)?;
    Ok(())
}

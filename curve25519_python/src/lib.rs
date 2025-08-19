use pyo3::prelude::*;
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::edwards::{CompressedEdwardsY};
use curve25519_dalek::constants::ED25519_BASEPOINT_POINT;

#[pyfunction]
fn scalar_exp(a: &[u8], x: u64) -> PyResult<Vec<u8>> {
    if a.len() != 32 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Scalar must be 32 bytes"));
    }
    let mut result = Scalar::ONE;
    let mut base = Scalar::from_bytes_mod_order(a.try_into().unwrap());
    let mut exp = x;
    while exp > 0 {
        if exp % 2 == 1 {
            result = result * base;
        }
        base = base * base;
        exp /= 2;
    }
    Ok(result.to_bytes().to_vec())
}
#[pyfunction]
fn scalar_multiply(scalar_bytes: &[u8]) -> PyResult<Vec<u8>> {
    if scalar_bytes.len() != 32 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Scalar must be 32 bytes"));
    }
    let scalar = Scalar::from_bytes_mod_order(scalar_bytes.try_into().unwrap());
    let point = scalar * ED25519_BASEPOINT_POINT;
    Ok(point.compress().to_bytes().to_vec())
}
#[pyfunction]
fn scalar_multiply_scalar(s1: &[u8], s2: &[u8]) -> PyResult<Vec<u8>> {
    if s1.len() != 32 || s2.len() != 32 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Scalars must be 32 bytes"));
    }
    let scalar1 = Scalar::from_bytes_mod_order(s1.try_into().unwrap());
    let scalar2 = Scalar::from_bytes_mod_order(s2.try_into().unwrap());
    let result = scalar1 * scalar2;
    Ok(result.to_bytes().to_vec())
}
#[pyfunction]
fn scalar_addition(s1: &[u8], s2: &[u8]) -> PyResult<Vec<u8>> {
    if s1.len() != 32 || s2.len() != 32 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Scalars must be 32 bytes"));
    }
    let scalar1 = Scalar::from_bytes_mod_order(s1.try_into().unwrap());
    let scalar2 = Scalar::from_bytes_mod_order(s2.try_into().unwrap());
    let result = scalar1 + scalar2;
    Ok(result.to_bytes().to_vec())
}
#[pyfunction]
fn scalar_subtraction(s1: &[u8], s2: &[u8]) -> PyResult<Vec<u8>> {
    if s1.len() != 32 || s2.len() != 32 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Scalars must be 32 bytes"));
    }
    let scalar1 = Scalar::from_bytes_mod_order(s1.try_into().unwrap());
    let scalar2 = Scalar::from_bytes_mod_order(s2.try_into().unwrap());
    let result = scalar1 - scalar2;
    Ok(result.to_bytes().to_vec())
}
#[pyfunction]
fn point_multiply(scalar_bytes: &[u8], point_bytes: &[u8]) -> PyResult<Vec<u8>> {
    if scalar_bytes.len() != 32 || point_bytes.len() != 32 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Inputs must be 32 bytes"));
    }
    let scalar = Scalar::from_bytes_mod_order(scalar_bytes.try_into().unwrap());
    let point = CompressedEdwardsY(point_bytes.try_into().unwrap())
        .decompress()
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid point"))?;
    let result_point = scalar * point;
    Ok(result_point.compress().to_bytes().to_vec())
}
#[pyfunction]
fn scalar_inverse(scalar_bytes: &[u8]) -> PyResult<Vec<u8>> {
    if scalar_bytes.len() != 32 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Scalar must be 32 bytes"));
    }
    let scalar = Scalar::from_bytes_mod_order(scalar_bytes.try_into().unwrap());
    if scalar == Scalar::ZERO {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Cannot invert zero"));
    }
    let inv = scalar.invert();
    Ok(inv.to_bytes().to_vec())
}
#[pyfunction]
fn point_addition(point1_bytes: &[u8], point2_bytes: &[u8]) -> PyResult<Vec<u8>> {
    if point1_bytes.len() != 32 || point2_bytes.len() != 32 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Points must be 32 bytes"));
    }
    let point1 = CompressedEdwardsY(point1_bytes.try_into().unwrap())
        .decompress()
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid point1"))?;
    let point2 = CompressedEdwardsY(point2_bytes.try_into().unwrap())
        .decompress()
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid point2"))?;
    let sum_point = point1 + point2;
    Ok(sum_point.compress().to_bytes().to_vec())
}
#[pyfunction]
fn point_subtraction(point1_bytes: &[u8], point2_bytes: &[u8]) -> PyResult<Vec<u8>> {
    if point1_bytes.len() != 32 || point2_bytes.len() != 32 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Points must be 32 bytes"));
    }
    let point1 = CompressedEdwardsY(point1_bytes.try_into().unwrap())
        .decompress()
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid point1"))?;
    let point2 = CompressedEdwardsY(point2_bytes.try_into().unwrap())
        .decompress()
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid point2"))?;
    let neg_point2 = -point2;
    let diff_point = point1 + neg_point2;
    Ok(diff_point.compress().to_bytes().to_vec())
}
#[pymodule]
fn curve25519_python<'py>(_py: Python<'py>, m: &Bound<'py, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(scalar_exp, m)?)?;
    m.add_function(wrap_pyfunction!(scalar_multiply, m)?)?;
    m.add_function(wrap_pyfunction!(scalar_multiply_scalar, m)?)?;
    m.add_function(wrap_pyfunction!(scalar_addition, m)?)?;
    m.add_function(wrap_pyfunction!(scalar_subtraction, m)?)?;
    m.add_function(wrap_pyfunction!(point_multiply, m)?)?;
    m.add_function(wrap_pyfunction!(scalar_inverse, m)?)?;
    m.add_function(wrap_pyfunction!(point_addition, m)?)?;
    m.add_function(wrap_pyfunction!(point_subtraction, m)?)?;
    Ok(())
}
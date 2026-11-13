//! Плотная линейная алгебра (на `nalgebra`): решение комплексных систем, матричная экспонента.

use nalgebra::DMatrix;
use num_complex::Complex64;

/// Dense complex matrix
pub type CMatrix = DMatrix<Complex64>;
/// Dense real matrix
pub type RMatrix = DMatrix<f64>;

/// Solution of the linear system `A X = B` with several right-hand sides (LU decomposition
/// with partial pivoting)
/// # Arguments
/// * `a` - square matrix `n × n`
/// * `b` - right-hand sides, `n × m`
/// # Result
/// Returns `Ok(X)`, `n × m`
/// # Errors
/// Returns descriptive string if the sizes do not match, the matrix is singular, or the
/// solution is not finite
pub fn solve(a: &CMatrix, b: &CMatrix) -> Result<CMatrix, String> {
    if a.nrows() != a.ncols() || a.nrows() != b.nrows() {
        return Err(format!(
            "linalg::solve(): matrix is {}x{}, right-hand side is {}x{}",
            a.nrows(),
            a.ncols(),
            b.nrows(),
            b.ncols()
        ));
    }
    let x = a
        .clone()
        .lu()
        .solve(b)
        .ok_or("linalg::solve(): the matrix is singular".to_string())?;
    if x.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()) {
        return Err("linalg::solve(): the matrix is singular (solution is not finite)".to_string());
    }
    Ok(x)
}

/// Matrix exponential `exp(A)` of a real square matrix (Padé approximation with scaling and
/// squaring)
/// # Errors
/// Returns descriptive string if the matrix is not square or the result is not finite
pub fn expm(a: &RMatrix) -> Result<RMatrix, String> {
    if a.nrows() != a.ncols() {
        return Err(format!(
            "linalg::expm(): matrix is {}x{}, must be square",
            a.nrows(),
            a.ncols()
        ));
    }
    if a.nrows() == 0 {
        return Ok(a.clone());
    }
    let result = a.exp();
    if result.iter().any(|v| !v.is_finite()) {
        return Err("linalg::expm(): the result is not finite".to_string());
    }
    Ok(result)
}

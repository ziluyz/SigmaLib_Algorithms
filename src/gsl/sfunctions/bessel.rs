mod gsl_signatures;

/// `GSL_EUNDRFLW`
const GSL_EUNDRFLW: i32 = 15;

/// Bessel function of the first kind `J_n(x)`
/// # Arguments
/// * `n` - integer order
/// * `x` - argument
/// # Result
/// Returns `Ok(J_n(x))`. Underflow is not an error: the result is the value returned by GSL
/// (zero).
/// # Errors
/// Returns descriptive string if the result is not finite (e.g. `x` is NaN) or GSL reports an
/// error (the latter works only after `gsl::errors::turn_off_error_handler()`, otherwise GSL
/// aborts the process)
pub fn jn(n: i32, x: f64) -> Result<f64, String> {
    let mut result = gsl_signatures::GslSfResult { val: 0.0, err: 0.0 };
    let status = unsafe { gsl_signatures::gsl_sf_bessel_Jn_e(n, x, &mut result) };
    if status == 0 || status == GSL_EUNDRFLW {
        // На нечисловом аргументе GSL возвращает «успех» и NaN.
        if result.val.is_finite() {
            Ok(result.val)
        } else {
            Err(format!(
                "Error in gsl_sf_bessel_Jn(n = {}, x = {}): result is not finite",
                n, x
            ))
        }
    } else {
        Err(format!(
            "Error in gsl_sf_bessel_Jn(n = {}, x = {}): {}",
            n,
            x,
            super::super::errors::get_error_description(status)
        ))
    }
}

/// Derivative of the Bessel function of the first kind `J_n'(x)`
/// (by the recurrence `J_n' = (J_{n-1} - J_{n+1}) / 2`)
/// # Arguments
/// * `n` - integer order
/// * `x` - argument
/// # Result
/// Returns `Ok(J_n'(x))`
/// # Errors
/// Returns the error of `jn`
pub fn djn(n: i32, x: f64) -> Result<f64, String> {
    Ok((jn(n - 1, x)? - jn(n + 1, x)?) * 0.5)
}

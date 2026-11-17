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

/// Largest `|x|` for which [`jn_djn`] sums the power series instead of calling GSL
pub const SERIES_LIMIT: f64 = 2.0;

/// `J_n(x)` and `J_n'(x)` together, fast for small arguments
///
/// For `|x| <= SERIES_LIMIT` the power series is summed to machine precision (no GSL call; the
/// terms decrease from the first one, there is no cancellation); for larger `|x|` the result is
/// that of [`jn`] and [`djn`]. Intended for inner loops such as the Larmor-radius factor
/// `J_n(k⊥ p⊥/H)` in the equations of electron motion.
/// # Arguments
/// * `n` - integer order
/// * `x` - argument
/// # Result
/// Returns `Ok((J_n(x), J_n'(x)))`
/// # Errors
/// Returns descriptive string if `x` is not finite, or the errors of [`jn`]
pub fn jn_djn(n: i32, x: f64) -> Result<(f64, f64), String> {
    if !x.is_finite() {
        return Err(format!("jn_djn(n = {n}, x = {x}): argument is not finite"));
    }
    if x.abs() > SERIES_LIMIT {
        return Ok((jn(n, x)?, djn(n, x)?));
    }
    // J_{-n} = (−1)^n J_n, J_n(−x) = (−1)^n J_n(x), J_n'(−x) = (−1)^{n+1} J_n'(x).
    let order = n.unsigned_abs();
    let sign_order = if n < 0 && order % 2 == 1 { -1.0 } else { 1.0 };
    let (j, dj) = series(order, x.abs());
    let odd = order % 2 == 1;
    let (j, dj) = if x < 0.0 {
        (if odd { -j } else { j }, if odd { dj } else { -dj })
    } else {
        (j, dj)
    };
    Ok((sign_order * j, sign_order * dj))
}

/// Series `J_n(x) = Σ_k (−1)^k (x/2)^{2k+n} / (k! (k+n)!)` and its derivative for `x >= 0`
fn series(n: u32, x: f64) -> (f64, f64) {
    if x == 0.0 {
        return (
            if n == 0 { 1.0 } else { 0.0 },
            if n == 1 { 0.5 } else { 0.0 },
        );
    }
    let half = 0.5 * x;
    let q = -half * half;
    // Первый член (x/2)^n / n!
    let mut term = 1.0;
    for i in 1..=n {
        term *= half / f64::from(i);
    }
    let (mut j, mut dj) = (0.0, 0.0);
    for k in 0..64u32 {
        j += term;
        // d/dx (x/2)^{2k+n} = (2k+n)/x · (x/2)^{2k+n}
        dj += term * f64::from(2 * k + n);
        let next = term * q / (f64::from(k + 1) * f64::from(k + 1 + n));
        if next.abs() <= f64::EPSILON * 0.25 * j.abs().max(f64::MIN_POSITIVE) {
            break;
        }
        term = next;
    }
    (j, dj / x)
}

/// `J_m(x)/x^m` — even entire function of `x`, regular at `x = 0` (value `1/(2^m m!)`)
///
/// Needed where a Bessel function of a small argument is divided by the matching power, e.g.
/// in the equations of electron motion written in the complex transverse momentum (the Larmor
/// radius may pass through zero). For `|x| <= SERIES_LIMIT` — power series to machine
/// precision, otherwise `jn(m, x)/x^m`.
/// # Arguments
/// * `m` - order, `m >= 0`
/// * `x` - argument
/// # Errors
/// Returns descriptive string if `x` is not finite, or the errors of [`jn`]
pub fn jn_over_power(m: u32, x: f64) -> Result<f64, String> {
    if !x.is_finite() {
        return Err(format!("jn_over_power(m = {m}, x = {x}): argument is not finite"));
    }
    let order = i32::try_from(m).map_err(|_| format!("jn_over_power(): order {m} is too large"))?;
    if x.abs() > SERIES_LIMIT {
        return Ok(jn(order, x)? / x.powi(order));
    }
    // Σ_k (−1)^k (x/2)^{2k} / (2^m k! (k+m)!)
    let mut term = 1.0;
    for i in 1..=m {
        term /= 2.0 * f64::from(i);
    }
    let q = -0.25 * x * x;
    let mut sum = 0.0;
    for k in 0..64u32 {
        sum += term;
        let next = term * q / (f64::from(k + 1) * f64::from(k + 1 + m));
        if next.abs() <= f64::EPSILON * 0.25 * sum.abs() {
            break;
        }
        term = next;
    }
    Ok(sum)
}

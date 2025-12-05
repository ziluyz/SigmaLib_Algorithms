/// Finding a root by bisection
/// # Arguments
/// * `func` - function of one argument
/// * `x_min` - left border
/// * `x_max` - right border
/// * `abs_err` - maximal absolute error
/// # Result
/// Returns `Ok(x)`, where `x` is a root of `func` on the interval `[x_min, x_max]` with accuracy
/// `abs_err`
/// # Errors
/// Returns descriptive string
pub fn find_root(
    func: &impl Fn(f64) -> Result<f64, String>,
    mut x_min: f64,
    mut x_max: f64,
    abs_err: f64,
) -> Result<f64, String> {
    if x_min >= x_max {
        return Err(format!(
            "xMin = {} > xMax = {} in find_root()",
            x_min, x_max
        ));
    }
    let mut f1 = func(x_min)?;
    if f1 * func(x_max)? > 0.0 {
        return Err(format!(
            "f(xMin = {}) and f(xMax = {}) have the same sign in find_root()",
            x_min, x_max
        ));
    }
    while x_max - x_min > abs_err {
        let x_mid = (x_min + x_max) / 2.0;
        let f_mid = func(x_mid)?;
        if f1 * f_mid < 0.0 {
            x_max = x_mid;
        } else {
            x_min = x_mid;
            f1 = f_mid;
        }
    }
    Ok((x_min + x_max) / 2.0)
}

/// Finding a sequence of roots
/// # Arguments
/// * `func` - function of one argument
/// * `x_start` - start point
/// * `min_root_dist` - minimal distance between roots
/// * `n` - number of roots
/// # Result
/// Returns `Ok(x)`, where `x` is a sequence of roots of `func` with maximal absolute error 1e-7
/// # Errors
/// Returns descriptive string
pub fn find_root_sequence(
    func: &impl Fn(f64) -> Result<f64, String>,
    x_start: f64,
    min_root_dist: f64,
    n: usize,
) -> Result<Vec<f64>, String> {
    let mut res = Vec::new();
    let mut x_min = x_start;
    let d = min_root_dist * 0.55;
    let x_step = min_root_dist * 0.5;
    for _ in 0..n {
        let mut a = x_min;
        let mut b = a + d;
        loop {
            match find_root(&func, a, b, 1e-7) {
                Ok(x) => {
                    res.push(x);
                    x_min = x + x_step;
                    break;
                }
                Err(_) => {
                    a += x_step;
                    b += x_step;
                }
            }
        }
    }
    Ok(res)
}

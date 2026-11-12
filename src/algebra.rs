pub mod argument_root;

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
    x_min: f64,
    x_max: f64,
    abs_err: f64,
) -> Result<f64, String> {
    if x_min >= x_max {
        return Err(format!(
            "xMin = {} > xMax = {} in find_root()",
            x_min, x_max
        ));
    }
    let f_min = func(x_min)?;
    if f_min == 0.0 {
        return Ok(x_min);
    }
    let f_max = func(x_max)?;
    if f_min * f_max > 0.0 {
        return Err(format!(
            "f(xMin = {}) and f(xMax = {}) have the same sign in find_root()",
            x_min, x_max
        ));
    }
    bisect(func, x_min, x_max, f_min, abs_err)
}

/// Бисекция на отрезке, на концах которого функция уже вычислена и имеет разные знаки
/// (`f1 = func(x_min)`, `f1 != 0`).
fn bisect(
    func: &impl Fn(f64) -> Result<f64, String>,
    mut x_min: f64,
    mut x_max: f64,
    mut f1: f64,
    abs_err: f64,
) -> Result<f64, String> {
    while x_max - x_min > abs_err {
        let x_mid = (x_min + x_max) / 2.0;
        let f_mid = func(x_mid)?;
        if f_mid == 0.0 {
            // Без этой проверки точное попадание в корень уводило поиск к x_max.
            return Ok(x_mid);
        }
        if f1 * f_mid < 0.0 {
            x_max = x_mid;
        } else {
            x_min = x_mid;
            f1 = f_mid;
        }
    }
    Ok((x_min + x_max) / 2.0)
}

/// Search limit of `find_root_sequence`: the maximal distance from the start point or from the
/// previous root to the next one, in units of `min_root_dist`.
pub const MAX_ROOT_GAP: f64 = 1e5;

/// Finding a sequence of roots
///
/// The search goes from `x_start` to the right in windows of `0.55 * min_root_dist` with the
/// step `0.5 * min_root_dist`; a root is located by bisection in a window where `func` changes
/// sign, and the search continues from `0.5 * min_root_dist` after it. Roots closer to each
/// other than `min_root_dist` can be missed.
/// # Arguments
/// * `func` - function of one argument
/// * `x_start` - start point
/// * `min_root_dist` - minimal distance between roots
/// * `n` - number of roots
/// # Result
/// Returns `Ok(x)`, where `x` is a sequence of `n` roots of `func` with maximal absolute error
/// 1e-7
/// # Errors
/// Returns descriptive string if `x_start` or `min_root_dist` is invalid, if the next root is
/// not found within `MAX_ROOT_GAP * min_root_dist` from the previous one (from `x_start` for
/// the first root), or the error of `func`
pub fn find_root_sequence(
    func: &impl Fn(f64) -> Result<f64, String>,
    x_start: f64,
    min_root_dist: f64,
    n: usize,
) -> Result<Vec<f64>, String> {
    let roots = scan_roots(func, x_start, f64::INFINITY, min_root_dist, n)?;
    if roots.len() < n {
        let from = roots.last().copied().unwrap_or(x_start);
        return Err(format!(
            "find_root_sequence(): found {} roots of {}: no sign change within {} after x = {}",
            roots.len(),
            n,
            MAX_ROOT_GAP * min_root_dist,
            from
        ));
    }
    Ok(roots)
}

/// Finding a sequence of roots on an interval
///
/// Same search as in `find_root_sequence`, but it stops at `x_max`, and running out of roots is
/// not an error.
/// # Arguments
/// * `func` - function of one argument
/// * `x_start` - start point
/// * `x_max` - right border of the search
/// * `min_root_dist` - minimal distance between roots
/// * `n` - maximal number of roots
/// # Result
/// Returns `Ok(x)`, where `x` is a sequence of the first roots of `func` on `[x_start, x_max]`
/// (at most `n`, possibly empty) with maximal absolute error 1e-7
/// # Errors
/// Returns descriptive string if the arguments are invalid, or the error of `func`
pub fn find_root_sequence_bounded(
    func: &impl Fn(f64) -> Result<f64, String>,
    x_start: f64,
    x_max: f64,
    min_root_dist: f64,
    n: usize,
) -> Result<Vec<f64>, String> {
    if !(x_max.is_finite() && x_max >= x_start) {
        return Err(format!(
            "xMax = {} must be finite and >= xStart = {} in find_root_sequence_bounded()",
            x_max, x_start
        ));
    }
    scan_roots(func, x_start, x_max, min_root_dist, n)
}

/// Общий поиск: до `n` корней правее `x_start`, не дальше `x_max` и не дальше
/// `MAX_ROOT_GAP * min_root_dist` от предыдущего корня.
fn scan_roots(
    func: &impl Fn(f64) -> Result<f64, String>,
    x_start: f64,
    x_max: f64,
    min_root_dist: f64,
    n: usize,
) -> Result<Vec<f64>, String> {
    if !x_start.is_finite() {
        return Err(format!(
            "xStart = {} is not finite in root sequence search",
            x_start
        ));
    }
    if !(min_root_dist.is_finite() && min_root_dist > 0.0) {
        return Err(format!(
            "minRootDist = {} must be finite and > 0 in root sequence search",
            min_root_dist
        ));
    }
    let mut res = Vec::new();
    let mut x_min = x_start;
    let d = min_root_dist * 0.55;
    let x_step = min_root_dist * 0.5;
    let max_windows = (MAX_ROOT_GAP / 0.5) as usize;
    'roots: for _ in 0..n {
        let mut a = x_min;
        let mut b = a + d;
        for _ in 0..max_windows {
            if a >= x_max {
                break 'roots;
            }
            let b_cut = b.min(x_max);
            let f_a = func(a)?;
            let x = if f_a == 0.0 {
                a
            } else {
                // Как и в find_root: одинаковые знаки — корня в окне нет.
                if f_a * func(b_cut)? > 0.0 {
                    a += x_step;
                    b += x_step;
                    continue;
                }
                bisect(func, a, b_cut, f_a, 1e-7)?
            };
            res.push(x);
            x_min = x + x_step;
            continue 'roots;
        }
        break;
    }
    Ok(res)
}

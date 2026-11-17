//! Newton's method for a map of the plane written as a function of a complex variable.

use super::argument_root::RectF64;
use num_complex::Complex64;

/// Settings of [`newton_2d`]
#[derive(Debug, Clone)]
pub struct Newton2dOptions {
    /// The point is accepted when `|f(z)| <= abs_tol`
    pub abs_tol: f64,
    /// The point is also accepted when the last step is not longer than `step_tol` in both
    /// coordinates (`re` — along the real axis, `im` — along the imaginary axis)
    pub step_tol: Complex64,
    /// Finite-difference steps of the Jacobian along both axes (positive; the sign is flipped
    /// near the upper boundary of the region)
    pub fd_step: Complex64,
    /// Maximal number of Newton iterations
    pub max_iterations: usize,
    /// Maximal number of step halvings per iteration
    pub max_halvings: usize,
}

/// A point found by [`newton_2d`]
#[derive(Debug, Clone, Copy)]
pub struct Newton2dRoot {
    pub z: Complex64,
    pub f: Complex64,
    /// Number of accepted Newton steps
    pub iterations: usize,
    /// Number of calls of the evaluator (each with 3 points)
    pub batches: usize,
}

/// Newton's method for a smooth map `R² → R²` written as `f: C → C` (`f` need not be
/// analytic: the full real 2×2 Jacobian is used)
///
/// The Jacobian is computed by forward differences. The evaluator receives three points at once
/// — `z`, `z + h_re`, `z + i·h_im` — and may compute them in parallel. A step that does not
/// decrease `|f|` or leaves `region` is halved. All points stay inside `region`.
/// # Arguments
/// * `eval` - batch evaluator: values of `f` at the given points, in the same order
/// * `start` - starting point inside `region`
/// * `region` - allowed region (trust region)
/// * `options` - tolerances and limits
/// # Result
/// Returns `Ok(Some(root))` when `|f| <= abs_tol` or the step became shorter than `step_tol`
/// (the caller decides whether `root.f` is small enough), `Ok(None)` if the iterations failed:
/// the limit of iterations, a degenerate Jacobian, no decrease of `|f|` after all halvings, or
/// a non-finite value
/// # Errors
/// Returns descriptive string if the arguments are invalid, the evaluator returns a wrong number
/// of values, or the error of the evaluator
pub fn newton_2d(
    eval: &impl Fn(&[Complex64]) -> Result<Vec<Complex64>, String>,
    start: Complex64,
    region: RectF64,
    options: &Newton2dOptions,
) -> Result<Option<Newton2dRoot>, String> {
    let h = options.fd_step;
    if !(h.re > 0.0 && h.im > 0.0 && h.re.is_finite() && h.im.is_finite()) {
        return Err(format!("newton_2d(): fd_step = {h} must be positive"));
    }
    if !(h.re < region.width() && h.im < region.height()) {
        return Err(format!(
            "newton_2d(): fd_step = {h} does not fit into the region {region:?}"
        ));
    }
    if !(options.abs_tol >= 0.0 && options.step_tol.re >= 0.0 && options.step_tol.im >= 0.0) {
        return Err("newton_2d(): tolerances must be non-negative".into());
    }
    if !inside(&region, start) {
        return Err(format!(
            "newton_2d(): start {start} is outside the region {region:?}"
        ));
    }
    let mut batches = 0;
    // Значение и якобиан (столбцы — производные по Re z и по Im z) в точке.
    let mut probe = |z: Complex64| -> Result<Option<(Complex64, Complex64, Complex64)>, String> {
        let hr = if z.re + h.re <= region.re_max { h.re } else { -h.re };
        let hi = if z.im + h.im <= region.im_max { h.im } else { -h.im };
        let points = [z, z + hr, z + Complex64::new(0.0, hi)];
        let values = eval(&points)?;
        batches += 1;
        if values.len() != 3 {
            return Err(format!(
                "newton_2d(): the evaluator returned {} values for 3 points",
                values.len()
            ));
        }
        if values.iter().any(|v| !v.is_finite()) {
            return Ok(None);
        }
        Ok(Some((
            values[0],
            (values[1] - values[0]) / hr,
            (values[2] - values[0]) / hi,
        )))
    };

    let Some((mut f, mut dfx, mut dfy)) = probe(start)? else {
        return Ok(None);
    };
    let mut z = start;
    for iteration in 0..=options.max_iterations {
        if f.norm() <= options.abs_tol {
            return Ok(Some(Newton2dRoot { z, f, iterations: iteration, batches }));
        }
        if iteration == options.max_iterations {
            break;
        }
        // J·d = −f: [Re dfx, Re dfy; Im dfx, Im dfy]·(d_re, d_im) = −(Re f, Im f).
        let det = dfx.re * dfy.im - dfy.re * dfx.im;
        if !(det.is_finite() && det != 0.0) {
            return Ok(None);
        }
        let d_re = (-f.re * dfy.im + f.im * dfy.re) / det;
        let d_im = (-dfx.re * f.im + dfx.im * f.re) / det;
        let mut step = Complex64::new(d_re, d_im);
        if !step.is_finite() {
            return Ok(None);
        }
        // Шаг уже в пределах допуска: дальше |f| убывать не обязан (ошибки округления).
        if step.re.abs() <= options.step_tol.re && step.im.abs() <= options.step_tol.im {
            return Ok(Some(Newton2dRoot { z, f, iterations: iteration, batches }));
        }
        let mut accepted = None;
        for _ in 0..=options.max_halvings {
            let candidate = z + step;
            if inside(&region, candidate)
                && let Some(probed) = probe(candidate)?
                && probed.0.norm() < f.norm()
            {
                accepted = Some((candidate, probed));
                break;
            }
            step *= 0.5;
        }
        let Some((candidate, probed)) = accepted else {
            return Ok(None);
        };
        z = candidate;
        (f, dfx, dfy) = probed;
        if step.re.abs() <= options.step_tol.re && step.im.abs() <= options.step_tol.im {
            return Ok(Some(Newton2dRoot { z, f, iterations: iteration + 1, batches }));
        }
    }
    Ok(None)
}

fn inside(r: &RectF64, z: Complex64) -> bool {
    z.re >= r.re_min && z.re <= r.re_max && z.im >= r.im_min && z.im <= r.im_max
}

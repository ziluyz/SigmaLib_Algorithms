//! Продолжение корня по параметру.

use crate::algebra::argument_root::RectF64;
use num_complex::Complex64;

/// Where to look for the root at the next value of the parameter
#[derive(Debug, Clone, Copy)]
pub struct Prediction {
    /// Predicted root (inside `window`)
    pub guess: Complex64,
    /// Search window (inside the global bounds)
    pub window: RectF64,
    /// `true` if the window has been widened after a failure at this parameter value
    pub widened: bool,
}

/// Following a root `z(p)` of an equation with a real parameter `p` (port of `mut::Follower`
/// of DeltaProject, with the failures reported)
///
/// For every parameter value the caller's `solve` looks for the root given a prediction:
/// the first point is searched in the window `first`; the second one in a window of the same size
/// around the first root; then the root is extrapolated through the last two (linear) or three
/// (quadratic, Newton's divided differences — the parameter steps may be uneven) roots, and the
/// window is `±3·|predicted change|` around the prediction, but not smaller than `min_half`.
/// If `solve` fails, the window is doubled around the prediction until it covers the global
/// bounds in both directions.
/// # Arguments
/// * `params` - parameter values, in the order of following (distinct)
/// * `first` - window for the first point (clipped to `global`)
/// * `global` - bounds that no window crosses
/// * `min_half` - minimal half-sizes of a window along Re and Im (positive)
/// * `solve` - `solve(index, p, prediction)`: `Ok(Some(z))` — the root, `Ok(None)` — not found in
///   the window
/// * `process` - called with every root found, in order
/// # Result
/// Returns `Ok(())` when all points are found
/// # Errors
/// Returns descriptive string if the arguments are invalid, the root is not found at some
/// parameter value (the points before it have been processed), or the errors of `solve`/`process`
pub fn follow(
    params: &[f64],
    first: RectF64,
    global: RectF64,
    min_half: Complex64,
    mut solve: impl FnMut(usize, f64, &Prediction) -> Result<Option<Complex64>, String>,
    mut process: impl FnMut(usize, f64, Complex64) -> Result<(), String>,
) -> Result<(), String> {
    if !(min_half.re > 0.0 && min_half.im > 0.0) {
        return Err(format!("follow(): min_half = {min_half} must be positive"));
    }
    if let Some(bad) = params.iter().position(|p| !p.is_finite()) {
        return Err(format!("follow(): parameter {bad} is not finite"));
    }
    let first = clip(&first, &global).ok_or_else(|| {
        format!("follow(): the first window {first:?} lies outside the bounds {global:?}")
    })?;
    // Найденные точки (p, z) — последние три.
    let mut found: Vec<(f64, Complex64)> = Vec::with_capacity(3);
    for (index, &p) in params.iter().enumerate() {
        let (guess, mut half) = match found[..] {
            [] => (center(&first), half_sizes(&first)),
            [(_, z0)] => (z0, half_sizes(&first)),
            [(p1, z1), (p2, z2)] => {
                let d12 = (z2 - z1) / (p2 - p1);
                let guess = z2 + d12 * (p - p2);
                (guess, spread(guess - z2))
            }
            [(p1, z1), (p2, z2), (p3, z3)] => {
                let d12 = (z2 - z1) / (p2 - p1);
                let d23 = (z3 - z2) / (p3 - p2);
                let d123 = (d23 - d12) / (p3 - p1);
                let guess = z3 + (p - p3) * (d23 + d123 * (p - p2));
                (guess, spread(guess - z3))
            }
            _ => unreachable!("at most three points are kept"),
        };
        if !guess.is_finite() {
            return Err(format!(
                "follow(): the prediction at p = {p} is not finite (repeated parameter value?)"
            ));
        }
        half = Complex64::new(half.re.max(min_half.re), half.im.max(min_half.im));
        let mut widened = false;
        let z = loop {
            let window = clip(
                &RectF64 {
                    re_min: guess.re - half.re,
                    re_max: guess.re + half.re,
                    im_min: guess.im - half.im,
                    im_max: guess.im + half.im,
                },
                &global,
            )
            .ok_or_else(|| {
                format!("follow(): the prediction {guess} at p = {p} lies outside the bounds {global:?}")
            })?;
            let guess = Complex64::new(
                guess.re.clamp(window.re_min, window.re_max),
                guess.im.clamp(window.im_min, window.im_max),
            );
            let prediction = Prediction { guess, window, widened };
            if let Some(z) = solve(index, p, &prediction)? {
                break z;
            }
            let covers = window.re_min <= global.re_min
                && window.re_max >= global.re_max
                && window.im_min <= global.im_min
                && window.im_max >= global.im_max;
            if covers {
                return Err(format!(
                    "follow(): no root found at p = {p} (point {index}) in the bounds {global:?}"
                ));
            }
            half *= 2.0;
            widened = true;
        };
        process(index, p, z)?;
        if found.len() == 3 {
            found.remove(0);
        }
        found.push((p, z));
    }
    Ok(())
}

fn center(r: &RectF64) -> Complex64 {
    Complex64::new(0.5 * (r.re_min + r.re_max), 0.5 * (r.im_min + r.im_max))
}

fn half_sizes(r: &RectF64) -> Complex64 {
    Complex64::new(0.5 * r.width(), 0.5 * r.height())
}

fn spread(change: Complex64) -> Complex64 {
    Complex64::new(3.0 * change.re.abs(), 3.0 * change.im.abs())
}

/// Intersection of rectangles, `None` if it is degenerate
fn clip(r: &RectF64, bounds: &RectF64) -> Option<RectF64> {
    RectF64::new(
        r.re_min.max(bounds.re_min),
        r.re_max.min(bounds.re_max),
        r.im_min.max(bounds.im_min),
        r.im_max.min(bounds.im_max),
    )
    .ok()
}

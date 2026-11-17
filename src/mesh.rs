//! Сетки для кусочно-линейных профилей.

/// What happens between node `i` and node `i + 1`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interval {
    /// The profile is linear on an interval of nonzero length
    Smooth,
    /// Both nodes have the same abscissa, the profile jumps from one value to another
    Jump,
}

/// Piecewise-linear profile resampled with a limited step
#[derive(Debug, Clone, PartialEq)]
pub struct Mesh {
    /// Abscissas of the nodes, non-decreasing; equal neighbours mark a jump
    pub x: Vec<f64>,
    /// `y[channel][node]`
    pub y: Vec<Vec<f64>>,
    /// `intervals[i]` describes the gap between nodes `i` and `i + 1`
    pub intervals: Vec<Interval>,
}

impl Mesh {
    /// Number of nodes
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// `true` if the mesh has no nodes (never the case for a mesh built by `refine`)
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Slope of a channel on the interval after node `i` (`None` for a jump or a bad index)
    pub fn slope(&self, channel: usize, i: usize) -> Option<f64> {
        if *self.intervals.get(i)? == Interval::Jump {
            return None;
        }
        let y = self.y.get(channel)?;
        Some((y.get(i + 1)? - y.get(i)?) / (self.x.get(i + 1)? - self.x.get(i)?))
    }
}

/// Resampling of piecewise-linear profiles given on a common grid
///
/// Every segment of the grid is divided into equal steps not longer than `max_step`; the
/// original points stay nodes of the mesh. Two consecutive points with the same abscissa are a
/// jump: both stay in the mesh as neighbouring nodes with `Interval::Jump` between them (if no
/// channel changes there, the duplicate point is dropped).
/// # Arguments
/// * `xs` - abscissas of the profile points, non-decreasing
/// * `ys` - channels: ordinates at the points of `xs`
/// * `max_step` - maximal distance between the nodes
/// # Result
/// Returns `Ok(Mesh)`
/// # Errors
/// Returns descriptive string if there are less than two points, the channels differ in length
/// from `xs`, `xs` decreases, the values are not finite, `max_step` is not positive, or the
/// profile has zero length
pub fn refine(xs: &[f64], ys: &[&[f64]], max_step: f64) -> Result<Mesh, String> {
    if xs.len() < 2 {
        return Err(format!(
            "mesh::refine(): need at least 2 points, got {}",
            xs.len()
        ));
    }
    if !(max_step > 0.0) {
        return Err(format!(
            "mesh::refine(): max_step = {max_step} must be positive"
        ));
    }
    for (channel, y) in ys.iter().enumerate() {
        if y.len() != xs.len() {
            return Err(format!(
                "mesh::refine(): channel {channel} has {} points, the grid has {}",
                y.len(),
                xs.len()
            ));
        }
        if let Some(bad) = y.iter().position(|v| !v.is_finite()) {
            return Err(format!(
                "mesh::refine(): channel {channel}, point {bad} is not a finite number"
            ));
        }
    }
    if let Some(bad) = xs.iter().position(|v| !v.is_finite()) {
        return Err(format!(
            "mesh::refine(): x[{bad}] is not a finite number"
        ));
    }
    let (first, last) = (xs[0], xs[xs.len() - 1]);
    if !(last > first) {
        return Err(format!(
            "mesh::refine(): the profile must have positive length, got [{first}; {last}]"
        ));
    }

    let mut mesh = Mesh {
        x: vec![first],
        y: ys.iter().map(|y| vec![y[0]]).collect(),
        intervals: Vec::new(),
    };
    for i in 1..xs.len() {
        let (x0, x1) = (xs[i - 1], xs[i]);
        if x1 < x0 {
            return Err(format!(
                "mesh::refine(): x decreases at point {i}: {x1} after {x0}"
            ));
        }
        if x1 == x0 {
            if ys.iter().any(|y| y[i] != y[i - 1]) {
                mesh.x.push(x1);
                for (out, y) in mesh.y.iter_mut().zip(ys) {
                    out.push(y[i]);
                }
                mesh.intervals.push(Interval::Jump);
            }
            continue;
        }
        // Небольшой запас, чтобы отрезок длиной ровно в n шагов не делился на n + 1.
        let steps = ((x1 - x0) / max_step * (1.0 - 1e-9)).ceil().max(1.0);
        if steps > 1e9 {
            return Err(format!(
                "mesh::refine(): max_step = {max_step} is too small for the segment [{x0}; {x1}]"
            ));
        }
        let steps = steps as usize;
        for j in 1..=steps {
            let t = j as f64 / steps as f64;
            // Конец отрезка — точно исходная точка, без ошибки округления.
            mesh.x.push(if j == steps { x1 } else { x0 + (x1 - x0) * t });
            for (out, y) in mesh.y.iter_mut().zip(ys) {
                out.push(if j == steps {
                    y[i]
                } else {
                    y[i - 1] + (y[i] - y[i - 1]) * t
                });
            }
            mesh.intervals.push(Interval::Smooth);
        }
    }
    Ok(mesh)
}

/// One-sided limits of a piecewise-linear profile at `x`
///
/// The profile is given by points `(xs[i], ys[i])` with non-decreasing `xs`; two points with the
/// same abscissa are a jump. Away from a jump both limits are the linearly interpolated value.
/// # Arguments
/// * `xs` - abscissas, non-decreasing, at least one point
/// * `ys` - ordinates at `xs`
/// * `x` - point, `xs[0] <= x <= xs[last]`
/// # Result
/// Returns `Ok((left, right))`: the limits from the left and from the right (at the ends of the
/// profile both are the end value)
/// # Errors
/// Returns descriptive string if the lengths differ, the profile is empty, `x` is not finite or
/// lies outside the profile
pub fn limits(xs: &[f64], ys: &[f64], x: f64) -> Result<(f64, f64), String> {
    if xs.len() != ys.len() || xs.is_empty() {
        return Err(format!(
            "mesh::limits(): need equal non-zero lengths, got {} abscissas and {} ordinates",
            xs.len(),
            ys.len()
        ));
    }
    let (first, last) = (xs[0], xs[xs.len() - 1]);
    if !x.is_finite() || x < first || x > last {
        return Err(format!(
            "mesh::limits(): x = {x} is outside the profile [{first}; {last}]"
        ));
    }
    let left = if x == first {
        ys[0]
    } else {
        // Первая точка с абсциссой >= x; перед ней точка левее x.
        let j = xs.partition_point(|&v| v < x);
        if xs[j] == x {
            ys[j]
        } else {
            lerp(xs[j - 1], ys[j - 1], xs[j], ys[j], x)
        }
    };
    let right = if x == last {
        ys[ys.len() - 1]
    } else {
        // Первая точка правее x; перед ней точка с абсциссой <= x.
        let j = xs.partition_point(|&v| v <= x);
        if xs[j - 1] == x {
            ys[j - 1]
        } else {
            lerp(xs[j - 1], ys[j - 1], xs[j], ys[j], x)
        }
    };
    Ok((left, right))
}

fn lerp(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    y0 + (y1 - y0) * (x - x0) / (x1 - x0)
}

/// Resampling of piecewise-linear profiles given on different grids
///
/// The counterpart of [`refine`] for profiles with their own abscissas (e.g. a cavity radius and
/// a magnetic field profile): the nodes of all profiles inside `[x_min, x_max]` and the ends of
/// the range are united, every profile is interpolated there (a jump of any profile becomes a
/// pair of nodes with the same abscissa), then the union grid is refined by [`refine`].
/// `mesh.y[i]` is the profile `profiles[i]`.
/// # Arguments
/// * `profiles` - pairs `(xs, ys)`; every profile must cover `[x_min, x_max]`
/// * `x_min`, `x_max` - range of the mesh, `x_min < x_max`
/// * `max_step` - maximal distance between the nodes
/// # Result
/// Returns `Ok(Mesh)`
/// # Errors
/// Returns descriptive string if there are no profiles, a profile is invalid (lengths,
/// decreasing or non-finite values) or does not cover the range, or the error of [`refine`]
pub fn refine_union(
    profiles: &[(&[f64], &[f64])],
    x_min: f64,
    x_max: f64,
    max_step: f64,
) -> Result<Mesh, String> {
    if profiles.is_empty() {
        return Err("mesh::refine_union(): no profiles".into());
    }
    if !(x_min.is_finite() && x_max.is_finite() && x_min < x_max) {
        return Err(format!(
            "mesh::refine_union(): invalid range [{x_min}; {x_max}]"
        ));
    }
    let mut knots = vec![x_min, x_max];
    for (p, (xs, ys)) in profiles.iter().enumerate() {
        if xs.len() != ys.len() || xs.is_empty() {
            return Err(format!(
                "mesh::refine_union(): profile {p} has {} abscissas and {} ordinates",
                xs.len(),
                ys.len()
            ));
        }
        if xs.iter().chain(ys.iter()).any(|v| !v.is_finite()) {
            return Err(format!(
                "mesh::refine_union(): profile {p} contains a non-finite value"
            ));
        }
        if let Some(i) = xs.windows(2).position(|w| w[1] < w[0]) {
            return Err(format!(
                "mesh::refine_union(): abscissas of profile {p} decrease at point {}",
                i + 1
            ));
        }
        if xs[0] > x_min || xs[xs.len() - 1] < x_max {
            return Err(format!(
                "mesh::refine_union(): profile {p} [{}; {}] does not cover the range [{x_min}; {x_max}]",
                xs[0],
                xs[xs.len() - 1]
            ));
        }
        knots.extend(xs.iter().copied().filter(|&x| x > x_min && x < x_max));
    }
    knots.sort_by(f64::total_cmp);
    knots.dedup();

    let mut grid = Vec::with_capacity(knots.len() * 2);
    let mut values: Vec<Vec<f64>> = vec![Vec::with_capacity(knots.len() * 2); profiles.len()];
    for &x in &knots {
        let lims = profiles
            .iter()
            .map(|(xs, ys)| limits(xs, ys, x))
            .collect::<Result<Vec<_>, _>>()?;
        // Внутри диапазона скачок даёт пару узлов; на концах берётся значение внутри диапазона.
        let jump = x > x_min && x < x_max && lims.iter().any(|(l, r)| l != r);
        if jump {
            grid.push(x);
            for (v, (l, _)) in values.iter_mut().zip(&lims) {
                v.push(*l);
            }
        }
        grid.push(x);
        for (v, (l, r)) in values.iter_mut().zip(&lims) {
            v.push(if x == x_max { *l } else { *r });
        }
    }
    let channels: Vec<&[f64]> = values.iter().map(Vec::as_slice).collect();
    refine(&grid, &channels, max_step)
}

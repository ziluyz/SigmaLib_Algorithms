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

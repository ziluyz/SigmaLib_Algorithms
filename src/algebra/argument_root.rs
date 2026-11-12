use num_complex::Complex64;
use std::collections::{HashMap, VecDeque};

pub type C64 = Complex64;

/// User-facing rectangle in C: Re in [re_min, re_max], Im in [im_min, im_max]
#[derive(Debug, Clone, Copy)]
pub struct RectF64 {
    pub re_min: f64,
    pub re_max: f64,
    pub im_min: f64,
    pub im_max: f64,
}

impl RectF64 {
    pub fn new(re_min: f64, re_max: f64, im_min: f64, im_max: f64) -> Result<Self, String> {
        if !re_min.is_finite() || !re_max.is_finite() || !im_min.is_finite() || !im_max.is_finite()
        {
            return Err("RectF64: bounds must be finite".into());
        }
        if re_min >= re_max || im_min >= im_max {
            return Err("RectF64: invalid bounds (min must be < max)".into());
        }
        Ok(Self {
            re_min,
            re_max,
            im_min,
            im_max,
        })
    }

    #[inline]
    pub fn width(&self) -> f64 {
        self.re_max - self.re_min
    }
    #[inline]
    pub fn height(&self) -> f64 {
        self.im_max - self.im_min
    }
}

/// Controls / limits.
#[derive(Debug, Clone)]
pub struct Limits {
    /// Maximal depth of the dyadic grid; a rectangle at this depth is reported as is
    /// (`StopReason::LimitsReached`).
    /// `max_depth + log2(samples_per_edge0) + max_winding_refine` must not exceed 58.
    pub max_depth: u32,
    /// Maximal number of processed rectangles; exceeding it is an error.
    pub max_rects: usize,

    /// Initial samples per edge (must be power of two).
    pub samples_per_edge0: u32,
    /// Max times to double samples_per_edge for winding computation.
    pub max_winding_refine: usize,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            max_depth: 40,
            max_rects: 200_000,
            samples_per_edge0: 32,
            max_winding_refine: 8,
        }
    }
}

/// Why we stopped refining a region.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopReason {
    AbsTolSatisfied,
    RelTolSatisfied,
    /// `Limits::max_depth` reached.
    LimitsReached,
    /// The rectangle could not be subdivided: for every tried position of the dividing lines
    /// the winding numbers of the parts were unreliable or inconsistent (zeros too close to
    /// each other for the given `Limits`, or `f` is dominated by rounding noise, as near a
    /// multiple zero). The rectangle is reported as is, `multiplicity` is the number of zeros
    /// inside.
    SubdivisionFailed,
}

/// Root report.
#[derive(Debug, Clone)]
pub struct Root {
    pub z: C64,
    pub fz: C64,
    pub abs_f: f64,

    /// Final localization rectangle (as f64 for reporting).
    pub rect: RectF64,
    /// Relative accuracy metric: diagonal length of rect.
    pub rect_diag: f64,

    /// Winding number on the final rectangle (zeros count with multiplicity).
    pub multiplicity: i64,

    pub stop_reason: StopReason,
    pub depth: u32,
}

/* -------------------- Dyadic grid + integer rectangles -------------------- */

/// Global coordinate system for dyadic keys:
/// Re = re0 + (re1-re0) * re_nom / 2^pow
/// Im = im0 + (im1-im0) * im_nom / 2^pow
///
/// We keep (nom, pow) normalized: nom==0 or odd.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Dyadic {
    pub nom: u64,
    pub pow: u32,
}

impl Dyadic {
    pub fn new(mut nom: u64, mut pow: u32) -> Self {
        if nom == 0 {
            return Self { nom: 0, pow: 0 };
        }
        while pow > 0 && (nom & 1) == 0 {
            nom >>= 1;
            pow -= 1;
        }
        Self { nom, pow }
    }
}

/// A cache key for a point z on the global dyadic grid of the initial rectangle.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PointKey {
    pub re: Dyadic,
    pub im: Dyadic,
}

/// Global mapping from dyadic to f64 coordinates.
#[derive(Debug, Clone, Copy)]
struct GlobalGrid {
    re0: f64,
    re_len: f64,
    im0: f64,
    im_len: f64,
}

impl GlobalGrid {
    fn from_rect(r: &RectF64) -> Self {
        Self {
            re0: r.re_min,
            re_len: r.re_max - r.re_min,
            im0: r.im_min,
            im_len: r.im_max - r.im_min,
        }
    }

    fn coord(&self, d: Dyadic, x0: f64, len: f64) -> f64 {
        if d.nom == 0 {
            return x0;
        }
        let denom = (1u128 << d.pow) as f64;
        x0 + len * (d.nom as f64) / denom
    }

    fn point(&self, key: PointKey) -> C64 {
        let re = self.coord(key.re, self.re0, self.re_len);
        let im = self.coord(key.im, self.im0, self.im_len);
        C64::new(re, im)
    }

    fn rect_f64(&self, r: &DRect) -> RectF64 {
        let re_min = self.coord(r.re_min_dyadic(), self.re0, self.re_len);
        let re_max = self.coord(r.re_max_dyadic(), self.re0, self.re_len);
        let im_min = self.coord(r.im_min_dyadic(), self.im0, self.im_len);
        let im_max = self.coord(r.im_max_dyadic(), self.im0, self.im_len);
        RectF64 {
            re_min,
            re_max,
            im_min,
            im_max,
        }
    }
}

/// Integer rectangle on dyadic grid at depth `depth`.
///
/// It represents:
///   re in [re_min_idx/2^depth, re_max_idx/2^depth] of the global Re segment
///   im in [im_min_idx/2^depth, im_max_idx/2^depth] of the global Im segment
///
/// Invariants:
/// - re_min_idx < re_max_idx, im_min_idx < im_max_idx
/// - indices are within [0, 2^depth]
#[derive(Debug, Clone, Copy)]
struct DRect {
    depth: u32,
    re_min_idx: u64,
    re_max_idx: u64,
    im_min_idx: u64,
    im_max_idx: u64,
}

impl DRect {
    fn root() -> Self {
        Self {
            depth: 0,
            re_min_idx: 0,
            re_max_idx: 1,
            im_min_idx: 0,
            im_max_idx: 1,
        }
    }

    #[inline]
    fn width_idx(&self) -> u64 {
        self.re_max_idx - self.re_min_idx
    }
    #[inline]
    fn height_idx(&self) -> u64 {
        self.im_max_idx - self.im_min_idx
    }

    /// Diagonal length in f64 coordinates.
    fn diag(&self, grid: &GlobalGrid) -> f64 {
        let rf = grid.rect_f64(self);
        (rf.width()).hypot(rf.height())
    }

    fn center_key(&self) -> PointKey {
        // center = (min+max)/2; in dyadic terms: nom = (min+max), pow = depth+1
        let re_nom = self.re_min_idx + self.re_max_idx;
        let im_nom = self.im_min_idx + self.im_max_idx;
        PointKey {
            re: Dyadic::new(re_nom, self.depth + 1),
            im: Dyadic::new(im_nom, self.depth + 1),
        }
    }

    fn re_min_dyadic(&self) -> Dyadic {
        Dyadic::new(self.re_min_idx, self.depth)
    }
    fn re_max_dyadic(&self) -> Dyadic {
        Dyadic::new(self.re_max_idx, self.depth)
    }
    fn im_min_dyadic(&self) -> Dyadic {
        Dyadic::new(self.im_min_idx, self.depth)
    }
    fn im_max_dyadic(&self) -> Dyadic {
        Dyadic::new(self.im_max_idx, self.depth)
    }

    /// Split into 4 children by the lines `re = re_min + width * re_eighths / 8` and
    /// `im = im_min + height * im_eighths / 8` (`4` = the middle).
    /// The middle split adds one level of depth, any other one adds three.
    /// Order: 0 bottom-left, 1 bottom-right, 2 top-left, 3 top-right.
    fn split4(&self, re_eighths: u64, im_eighths: u64) -> [DRect; 4] {
        let (shift, re_part, im_part) = if re_eighths == 4 && im_eighths == 4 {
            (1, 1, 1)
        } else {
            (3, re_eighths, im_eighths)
        };
        let d1 = self.depth + shift;
        let re0 = self.re_min_idx << shift;
        let re1 = self.re_max_idx << shift;
        let im0 = self.im_min_idx << shift;
        let im1 = self.im_max_idx << shift;

        let re_mid = re0 + self.width_idx() * re_part;
        let im_mid = im0 + self.height_idx() * im_part;

        let child = |re_min_idx, re_max_idx, im_min_idx, im_max_idx| DRect {
            depth: d1,
            re_min_idx,
            re_max_idx,
            im_min_idx,
            im_max_idx,
        };
        [
            child(re0, re_mid, im0, im_mid),
            child(re_mid, re1, im0, im_mid),
            child(re0, re_mid, im_mid, im1),
            child(re_mid, re1, im_mid, im1),
        ]
    }
}

/* -------------------- Cache -------------------- */

#[derive(Debug, Default)]
pub struct EvalCache {
    map: HashMap<PointKey, C64>,
    pub hits: usize,
    pub misses: usize,
}

impl EvalCache {
    pub fn new() -> Self {
        Self::default()
    }

    fn get_or_eval<F>(&mut self, key: PointKey, z: C64, f: &F) -> C64
    where
        F: Fn(C64) -> C64,
    {
        if let Some(v) = self.map.get(&key) {
            self.hits += 1;
            return *v;
        }
        self.misses += 1;
        let v = f(z);
        self.map.insert(key, v);
        v
    }
}

/* -------------------- Main API -------------------- */

/// Positions of the dividing lines (in eighths of the rectangle, Re and Im) tried in turn when
/// a rectangle is subdivided. The middle goes first; the others move the lines away from zeros
/// that happen to lie on them or too close to them.
const SPLITS: [(u64, u64); 5] = [(4, 4), (3, 5), (5, 3), (3, 3), (5, 5)];

/// Doublings of the sample count allowed in the first, cheap round over `SPLITS`: moving the
/// dividing lines away from a zero is cheaper than resolving a contour passing next to it.
const CHEAP_WINDING_REFINE: usize = 2;

/// Finding zeros of an analytic function in a rectangle of the complex plane by the argument
/// principle with adaptive subdivision.
///
/// The number of zeros inside a rectangle is the winding number of `f` along its boundary.
/// A rectangle with a nonzero winding number is divided into four parts until one of the stop
/// criteria is met; then its center is reported as a zero.
/// # Arguments
/// * `f` - function, analytic in `rect` (a pole is counted as a zero of negative multiplicity)
/// * `rect` - search area; no zero may lie on its boundary
/// * `abs_tol` - a rectangle with a single zero is accepted when `|f(center)| <= abs_tol`
///   (0 turns this criterion off)
/// * `rel_tol` - a rectangle is accepted when its diagonal is `<= rel_tol`
/// * `stop_after_first` - return as soon as one zero is found
/// * `limits` - limits of the search
/// # Result
/// Returns `Ok(roots)`: one `Root` per located zero or cluster of zeros not separated at the
/// requested accuracy; `Root::multiplicity` is the number of zeros in `Root::rect`, counted with
/// multiplicity, `Root::stop_reason` tells why the rectangle was accepted. The sum of
/// multiplicities equals the number of zeros in `rect`.
/// # Errors
/// Returns descriptive string if arguments are invalid, if the winding number along the boundary
/// of `rect` cannot be computed reliably (a zero on the boundary or too close to it — move the
/// boundary), or if `limits.max_rects` is exceeded
pub fn find_zeros_by_argument<F>(
    f: F,
    rect: RectF64,
    abs_tol: f64,
    rel_tol: f64,
    stop_after_first: bool,
    limits: Limits,
) -> Result<Vec<Root>, String>
where
    F: Fn(C64) -> C64,
{
    if !(abs_tol.is_finite() && abs_tol >= 0.0) {
        return Err("abs_tol must be finite and >= 0".into());
    }
    if !(rel_tol.is_finite() && rel_tol > 0.0) {
        return Err("rel_tol must be finite and > 0".into());
    }
    if limits.samples_per_edge0 < 2 || !limits.samples_per_edge0.is_power_of_two() {
        return Err("limits.samples_per_edge0 must be a power of two and >= 2".into());
    }
    // Индексы точек на сетке — u64: глубина вместе с дроблением контура должна в него помещаться
    // (с запасом на сдвинутое деление, которое добавляет три уровня глубины сразу).
    let max_bits = limits.max_depth as u64
        + limits.samples_per_edge0.trailing_zeros() as u64
        + limits.max_winding_refine as u64;
    if max_bits > 58 {
        return Err(
            "limits: max_depth + log2(samples_per_edge0) + max_winding_refine must be <= 58".into(),
        );
    }

    let grid = GlobalGrid::from_rect(&rect);
    let mut cache = EvalCache::new();

    let mut roots: Vec<Root> = Vec::new();

    // Start from integer root rectangle.
    let r0 = DRect::root();

    let n0 = match winding_number_rect(
        &f,
        &r0,
        limits.max_winding_refine,
        &limits,
        &grid,
        &mut cache,
    ) {
        Some(n) => n,
        None => {
            return Err(
                "Cannot compute the winding number on the boundary of the rectangle: a zero \
                 (or a singularity) lies on the boundary or too close to it; move the boundary"
                    .into(),
            );
        }
    };

    if n0 == 0 {
        return Ok(roots);
    }

    #[derive(Debug, Clone)]
    struct Cand {
        r: DRect,
        n: i64,
    }

    let mut queue = VecDeque::new();
    queue.push_back(Cand { r: r0, n: n0 });

    let mut processed = 0usize;

    while let Some(c) = queue.pop_front() {
        processed += 1;
        if processed > limits.max_rects {
            return Err("Rectangle processing budget exceeded".into());
        }

        let center_key = c.r.center_key();
        let center_z = grid.point(center_key);
        let f_center = cache.get_or_eval(center_key, center_z, &f); // cache center too
        let abs_f = f_center.norm();
        let diag = c.r.diag(&grid);

        // Stop criteria: OR (less strict triggers first).
        // Малое |f| в центре говорит о найденном нуле, только если ноль в прямоугольнике один:
        // иначе остальные были бы потеряны.
        let mut stop_reason = if c.n.abs() == 1 && abs_f <= abs_tol {
            Some(StopReason::AbsTolSatisfied)
        } else if diag <= rel_tol {
            Some(StopReason::RelTolSatisfied)
        } else if c.r.depth >= limits.max_depth {
            Some(StopReason::LimitsReached)
        } else {
            None
        };

        if stop_reason.is_none() {
            match subdivide(&f, &c.r, c.n, &limits, &grid, &mut cache) {
                Some(children) => {
                    for (r, n) in children {
                        queue.push_back(Cand { r, n });
                    }
                    continue;
                }
                None => stop_reason = Some(StopReason::SubdivisionFailed),
            }
        }

        if let Some(reason) = stop_reason {
            roots.push(Root {
                z: center_z,
                fz: f_center,
                abs_f,
                rect: grid.rect_f64(&c.r),
                rect_diag: diag,
                multiplicity: c.n,
                stop_reason: reason,
                depth: c.r.depth,
            });

            if stop_after_first {
                return Ok(roots);
            }
        }
    }

    Ok(roots)
}

/// Divides `r` (winding number `n`) into four parts and returns those with a nonzero winding
/// number. A division is accepted if the winding numbers of all four parts are reliable and sum
/// up to `n`. Returns `None` if no position of the dividing lines from `SPLITS` is accepted.
fn subdivide<F>(
    f: &F,
    r: &DRect,
    n: i64,
    limits: &Limits,
    grid: &GlobalGrid,
    cache: &mut EvalCache,
) -> Option<Vec<(DRect, i64)>>
where
    F: Fn(C64) -> C64,
{
    let cheap = limits.max_winding_refine.min(CHEAP_WINDING_REFINE);
    let budgets = if cheap < limits.max_winding_refine {
        vec![cheap, limits.max_winding_refine]
    } else {
        vec![cheap]
    };

    for max_refine in budgets {
        'split: for (re_eighths, im_eighths) in SPLITS {
            let mut children = Vec::with_capacity(4);
            let mut total = 0;
            for child in r.split4(re_eighths, im_eighths) {
                let Some(wn) = winding_number_rect(f, &child, max_refine, limits, grid, cache)
                else {
                    continue 'split;
                };
                total += wn;
                if wn != 0 {
                    children.push((child, wn));
                }
            }
            if total == n {
                return Some(children);
            }
        }
    }
    None
}

/* -------------------- Winding number on integer rectangles -------------------- */

/// Winding number of `f` along the boundary of `r`, or `None` if it cannot be computed reliably
/// with at most `max_refine` doublings of the number of samples.
fn winding_number_rect<F>(
    f: &F,
    r: &DRect,
    max_refine: usize,
    limits: &Limits,
    grid: &GlobalGrid,
    cache: &mut EvalCache,
) -> Option<i64>
where
    F: Fn(C64) -> C64,
{
    let mut samples = limits.samples_per_edge0;

    for _ in 0..=max_refine {
        match winding_pass(f, r, samples, grid, cache) {
            WindingPass::Smooth(n) => return Some(n),
            WindingPass::Rough => {}
            // Точки более частой выборки включают текущие: дробить бесполезно.
            WindingPass::Singular => return None,
        }
        samples *= 2;
    }

    None
}

enum WindingPass {
    /// `log f` changes slowly between adjacent samples: the winding number is reliable.
    Smooth(i64),
    /// Sampling is too coarse.
    Rough,
    /// `f` is zero or not finite at a sample point.
    Singular,
}

/// Maximal change of `arg f` and of `ln |f|` between adjacent samples of a reliable contour.
///
/// Обе величины — части приращения `log f`. Одного контроля фазы мало: рядом с m-кратным нулём,
/// лежащим почти на контуре, фаза между соседними точками меняется почти на `m·π`, и при чётном
/// `m` это неотличимо от малого изменения. Модуль при этом меняется между соседними точками не
/// менее чем в `3^m` раз, что и ловит вторая проверка.
const LOG_STEP_MAX: f64 = 1.0;

fn winding_pass<F>(
    f: &F,
    r: &DRect,
    samples_per_edge: u32,
    grid: &GlobalGrid,
    cache: &mut EvalCache,
) -> WindingPass
where
    F: Fn(C64) -> C64,
{
    // Unique contour sampling: n points per edge, the end corner of an edge belongs to the next
    // one.
    let n = samples_per_edge;

    let mut w: Vec<C64> = Vec::with_capacity(4 * n as usize);
    for edge in 0..4 {
        for k in 0..n {
            let key = boundary_key(r, edge, k, n);
            let fz = cache.get_or_eval(key, grid.point(key), f);
            if !fz.re.is_finite() || !fz.im.is_finite() || (fz.re == 0.0 && fz.im == 0.0) {
                return WindingPass::Singular;
            }
            w.push(fz);
        }
    }

    // Winding number via phase increments.
    let m = w.len();
    let mut total = 0.0f64;
    let mut smooth = true;

    for i in 0..m {
        let a = w[i];
        let b = w[(i + 1) % m];
        let prod = a.conj() * b;
        let delta = prod.im.atan2(prod.re);
        total += delta;
        let log_ratio = (b.norm() / a.norm()).ln();
        // Сравнение в такой форме отбраковывает и NaN.
        if !(delta.abs() <= LOG_STEP_MAX && log_ratio.abs() <= LOG_STEP_MAX) {
            smooth = false;
        }
    }

    if smooth {
        WindingPass::Smooth((total / (2.0 * std::f64::consts::PI)).round() as i64)
    } else {
        WindingPass::Rough
    }
}

/// Boundary point key for integer rectangle `r` at depth `d`.
///
/// We represent any boundary point as dyadics in the *global* coordinates:
/// - rectangle bounds are at denominator 2^d
/// - sampling adds denominator 2^s where samples_per_edge=2^s
/// => varying coordinate has denominator 2^(d+s)
///
/// For each edge:
/// 0 bottom: im fixed at im_min, re varies min->max
/// 1 right : re fixed at re_max, im varies min->max
/// 2 top   : im fixed at im_max, re varies max->min
/// 3 left  : re fixed at re_min, im varies max->min
fn boundary_key(r: &DRect, edge: u8, k: u32, samples_per_edge: u32) -> PointKey {
    let d = r.depth;
    let s = samples_per_edge.trailing_zeros(); // since power of two
    let pow_var = d + s;

    // Convert rectangle indices to the same denominator 2^(d+s) by shifting left by s.
    let re0 = r.re_min_idx << s;
    let re1 = r.re_max_idx << s;
    let im0 = r.im_min_idx << s;
    let im1 = r.im_max_idx << s;

    // Parameter t = k / samples_per_edge, k in [0..samples_per_edge-1].
    // idx = start + (end-start)*k/samples_per_edge; the bounds are already multiplied by
    // samples_per_edge, so the step is the width (height) of the rectangle in its own indices.
    let kk = k as u64;
    let re_step = r.width_idx();
    let im_step = r.height_idx();

    let (re_nom, im_nom) = match edge {
        0 => (re0 + re_step * kk, im0),
        1 => (re1, im0 + im_step * kk),
        2 => (re1 - re_step * kk, im1),
        3 => (re0, im1 - im_step * kk),
        _ => unreachable!(),
    };

    PointKey {
        re: Dyadic::new(re_nom, pow_var),
        im: Dyadic::new(im_nom, pow_var),
    }
}

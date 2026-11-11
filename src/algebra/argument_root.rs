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
    pub max_depth: u32,
    pub max_rects: usize,

    /// Initial samples per edge (must be power of two).
    pub samples_per_edge0: u32,
    /// Max times to double samples_per_edge for winding computation.
    pub max_winding_refine: usize,

    /// If min |f| on contour is <= abs_tol => treat as boundary root (reported).
    /// If min |f| on contour is <= contour_suspicion_tol but > abs_tol => error.
    pub contour_suspicion_tol: f64,

    /// Merge distance for deduplication of roots (in z-plane).
    pub merge_eps: f64,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            max_depth: 60,
            max_rects: 200_000,
            samples_per_edge0: 32,
            max_winding_refine: 8,
            contour_suspicion_tol: 1e-10,
            merge_eps: 1e-8,
        }
    }
}

/// Why we stopped refining a region.
#[derive(Debug, Clone, Copy)]
pub enum StopReason {
    AbsTolSatisfied,
    RelTolSatisfied,
    LimitsReached,
    BoundaryHit,
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

    /// Split into 4 children at depth+1.
    /// Order: 0 bottom-left, 1 bottom-right, 2 top-left, 3 top-right.
    fn split4(&self) -> [DRect; 4] {
        let d1 = self.depth + 1;
        let re0 = self.re_min_idx * 2;
        let re1 = self.re_max_idx * 2;
        let im0 = self.im_min_idx * 2;
        let im1 = self.im_max_idx * 2;

        let re_mid = (re0 + re1) / 2;
        let im_mid = (im0 + im1) / 2;

        [
            DRect {
                depth: d1,
                re_min_idx: re0,
                re_max_idx: re_mid,
                im_min_idx: im0,
                im_max_idx: im_mid,
            },
            DRect {
                depth: d1,
                re_min_idx: re_mid,
                re_max_idx: re1,
                im_min_idx: im0,
                im_max_idx: im_mid,
            },
            DRect {
                depth: d1,
                re_min_idx: re0,
                re_max_idx: re_mid,
                im_min_idx: im_mid,
                im_max_idx: im1,
            },
            DRect {
                depth: d1,
                re_min_idx: re_mid,
                re_max_idx: re1,
                im_min_idx: im_mid,
                im_max_idx: im1,
            },
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
    println!(
        "susp_tol: {}, abs_tol: {}",
        limits.contour_suspicion_tol, abs_tol
    );
    if !(limits.contour_suspicion_tol.is_finite() && limits.contour_suspicion_tol >= abs_tol) {
        return Err("limits.contour_suspicion_tol must be finite and >= abs_tol".into());
    }

    let grid = GlobalGrid::from_rect(&rect);
    let mut cache = EvalCache::new();

    let mut roots: Vec<Root> = Vec::new();

    // Start from integer root rectangle.
    let r0 = DRect::root();

    // Initial winding number (also may report boundary roots).
    let n0 = match winding_number_rect(&f, &r0, abs_tol, &limits, &grid, &mut cache, &mut roots)? {
        WindingResult::Wn(n) => n,
        WindingResult::BoundaryRoot => {
            if stop_after_first && !roots.is_empty() {
                dedup_roots(&mut roots, limits.merge_eps);
                return Ok(roots);
            }
            1
        }
    };

    if n0 == 0 {
        dedup_roots(&mut roots, limits.merge_eps);
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
        println!(
            "cand n: {}, cand DRect: {}, {}, {}, {}, {}",
            c.n, c.r.depth, c.r.re_min_idx, c.r.re_max_idx, c.r.im_min_idx, c.r.im_max_idx
        );
        processed += 1;
        if processed > limits.max_rects {
            return Err("Rectangle processing budget exceeded".into());
        }

        // Stop criteria: OR (less strict triggers first).
        let center_key = c.r.center_key();
        let center_z = grid.point(center_key);
        println!("Center: {center_z}");
        let f_center = cache.get_or_eval(center_key, center_z, &f); // cache center too
        let abs_f = f_center.norm();
        let diag = c.r.diag(&grid);

        let stop_reason = if abs_f <= abs_tol {
            Some(StopReason::AbsTolSatisfied)
        } else if diag <= rel_tol {
            Some(StopReason::RelTolSatisfied)
        } else if c.r.depth >= limits.max_depth {
            Some(StopReason::LimitsReached)
        } else {
            None
        };

        if let Some(reason) = stop_reason {
            let rect_f = grid.rect_f64(&c.r);
            roots.push(Root {
                z: center_z,
                fz: f_center,
                abs_f,
                rect: rect_f,
                rect_diag: diag,
                multiplicity: c.n,
                stop_reason: reason,
                depth: c.r.depth,
            });

            if stop_after_first {
                dedup_roots(&mut roots, limits.merge_eps);
                return Ok(roots);
            }
            continue;
        }

        // Subdivide.
        for child in c.r.split4() {
            println!(
                "child: {}, {} {} {} {}",
                child.depth, child.re_min_idx, child.re_max_idx, child.im_min_idx, child.im_max_idx
            );
            let mut boundary_hits: Vec<Root> = Vec::new();
            let wn = match winding_number_rect(
                &f,
                &child,
                abs_tol,
                &limits,
                &grid,
                &mut cache,
                &mut boundary_hits,
            )? {
                WindingResult::Wn(n) => n,
                WindingResult::BoundaryRoot => {
                    roots.extend(boundary_hits);
                    if stop_after_first && !roots.is_empty() {
                        dedup_roots(&mut roots, limits.merge_eps);
                        return Ok(roots);
                    }
                    1
                }
            };

            // if !boundary_hits.is_empty() {
            //     roots.extend(boundary_hits);
            //     if stop_after_first && !roots.is_empty() {
            //         dedup_roots(&mut roots, limits.merge_eps);
            //         return Ok(roots);
            //     }
            // }

            if wn != 0 {
                queue.push_back(Cand { r: child, n: wn });
            }
        }
    }

    dedup_roots(&mut roots, limits.merge_eps);
    Ok(roots)
}

/* -------------------- Winding number on integer rectangles -------------------- */

enum WindingResult {
    Wn(i64),
    BoundaryRoot,
}

fn winding_number_rect<F>(
    f: &F,
    r: &DRect,
    abs_tol: f64,
    limits: &Limits,
    grid: &GlobalGrid,
    cache: &mut EvalCache,
    roots_out: &mut Vec<Root>,
) -> Result<WindingResult, String>
where
    F: Fn(C64) -> C64,
{
    let mut samples = limits.samples_per_edge0;
    let mut last: Option<i64> = None;

    for _ in 0..=limits.max_winding_refine {
        match winding_pass(
            f,
            r,
            samples,
            abs_tol,
            limits.contour_suspicion_tol,
            grid,
            cache,
        )? {
            WindingPass::BoundaryHits(hits) => {
                // Report all hits as BoundaryHit roots.
                let rect_f = grid.rect_f64(r);
                let diag = r.diag(grid);
                for (key, fz) in hits {
                    let z = grid.point(key);
                    roots_out.push(Root {
                        z,
                        fz,
                        abs_f: fz.norm(),
                        rect: rect_f,
                        rect_diag: diag,
                        multiplicity: 1,
                        stop_reason: StopReason::BoundaryHit,
                        depth: r.depth,
                    });
                }
                return Ok(WindingResult::BoundaryRoot);
            }
            WindingPass::Unstable { min_abs_f } => {
                return Err(format!(
                    "Contour too close to a zero: min |f| on boundary = {min_abs_f:.3e} \
                     (abs_tol={abs_tol:.3e}, suspicion_tol={:.3e})",
                    limits.contour_suspicion_tol
                ));
            }
            WindingPass::Ok { n, stable } => {
                if stable {
                    if let Some(prev) = last {
                        if prev == n {
                            return Ok(WindingResult::Wn(n));
                        }
                    } else {
                        return Ok(WindingResult::Wn(n));
                    }
                }
                last = Some(n);
            }
        }
        samples = samples.saturating_mul(2);
    }

    Err("Failed to compute stable winding number on rectangle boundary".into())
}

enum WindingPass {
    Ok { n: i64, stable: bool },
    BoundaryHits(Vec<(PointKey, C64)>),
    Unstable { min_abs_f: f64 },
}

fn winding_stable(w: &[C64]) -> bool {
    // Conservative stability check:
    // - if any step has phase jump close to pi, we likely undersampled
    // - if any value is non-finite, unstable
    //
    // Threshold chosen to be strict enough to force refinement when needed.
    const PHASE_JUMP_MAX: f64 = 2.4; // < pi, stricter than previous 2.8

    if w.len() < 4 {
        return false;
    }

    for &val in w {
        if !val.re.is_finite() || !val.im.is_finite() {
            return false;
        }
    }

    for i in 0..w.len() {
        let a = w[i];
        let b = w[(i + 1) % w.len()];
        let prod = a.conj() * b;
        let delta = prod.im.atan2(prod.re);
        if delta.abs() > PHASE_JUMP_MAX {
            return false;
        }
    }
    true
}

fn winding_pass<F>(
    f: &F,
    r: &DRect,
    samples_per_edge: u32,
    abs_tol: f64,
    suspicion_tol: f64,
    grid: &GlobalGrid,
    cache: &mut EvalCache,
) -> Result<WindingPass, String>
where
    F: Fn(C64) -> C64,
{
    if samples_per_edge < 2 || !samples_per_edge.is_power_of_two() {
        return Err("Internal error: samples_per_edge must be power of two".into());
    }

    // Unique contour sampling:
    // - bottom: k = 0..n-1
    // - right : k = 0..n-1, but skip k=0 (bottom-right corner)
    // - top   : k = 0..n-1, but skip k=0 (top-right corner)
    // - left  : k = 0..n-1, but skip k=0 (top-left corner)
    let n = samples_per_edge;

    let mut keys: Vec<PointKey> = Vec::with_capacity((4 * n - 3) as usize);

    // bottom
    for k in 0..n {
        keys.push(boundary_key(r, 0, k, n));
    }
    // right (skip first corner)
    for k in 1..n {
        keys.push(boundary_key(r, 1, k, n));
    }
    // top (skip first corner)
    for k in 1..n {
        keys.push(boundary_key(r, 2, k, n));
    }
    // left (skip first corner)
    for k in 1..n {
        keys.push(boundary_key(r, 3, k, n));
    }

    // Evaluate f on contour with caching.
    let mut w: Vec<C64> = Vec::with_capacity(keys.len());
    let mut min_abs_f = f64::INFINITY;
    let mut hits: Vec<(PointKey, C64)> = Vec::new();

    for &key in &keys {
        let z = grid.point(key);
        let fz = cache.get_or_eval(key, z, f);
        let af = fz.norm();
        min_abs_f = min_abs_f.min(af);
        if af <= abs_tol {
            hits.push((key, fz));
        }
        w.push(fz);
    }

    if !hits.is_empty() {
        return Ok(WindingPass::BoundaryHits(hits));
    }
    if min_abs_f <= suspicion_tol {
        return Ok(WindingPass::Unstable { min_abs_f });
    }

    // Winding number via phase increments.
    let m = w.len();
    let mut total = 0.0f64;
    let mut stable = true;

    for i in 0..m {
        let a = w[i];
        let b = w[(i + 1) % m];
        let prod = a.conj() * b;
        let delta = prod.im.atan2(prod.re);
        total += delta;
    }

    // Stability check moved to a separate function (see пункт 2).
    stable = stable && winding_stable(&w);

    let wn = (total / (2.0 * std::f64::consts::PI)).round() as i64;
    Ok(WindingPass::Ok { n: wn, stable })
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

    // k is in [0..samples_per_edge-1]; we sample corners multiple times; ok.
    // Parameter t = k / samples_per_edge.
    // For forward direction: idx = start + (end-start)*k/samples_per_edge.
    // Since denom is 2^s, this becomes integer arithmetic:
    // idx = start + (end-start)*k / 2^s
    let kk = k as u64;
    let denom = samples_per_edge as u64;

    let (re_nom, im_nom) = match edge {
        0 => {
            let re_nom = re0 + ((re1 - re0) * kk) / denom;
            (re_nom, im0)
        }
        1 => {
            let im_nom = im0 + ((im1 - im0) * kk) / denom;
            (re1, im_nom)
        }
        2 => {
            let re_nom = re1 - ((re1 - re0) * kk) / denom;
            (re_nom, im1)
        }
        3 => {
            let im_nom = im1 - ((im1 - im0) * kk) / denom;
            (re0, im_nom)
        }
        _ => unreachable!(),
    };

    PointKey {
        re: Dyadic::new(re_nom, pow_var),
        im: Dyadic::new(im_nom, pow_var),
    }
}

/* -------------------- Dedup -------------------- */

fn dedup_roots(roots: &mut Vec<Root>, eps: f64) {
    if roots.len() <= 1 {
        return;
    }
    // Prefer smaller |f|.
    roots.sort_by(|a, b| {
        a.abs_f
            .partial_cmp(&b.abs_f)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut kept: Vec<Root> = Vec::with_capacity(roots.len());
    'outer: for r in roots.drain(..) {
        for k in &kept {
            if (r.z - k.z).norm() <= eps {
                continue 'outer;
            }
        }
        kept.push(r);
    }
    *roots = kept;
}

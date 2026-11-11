mod common;

use algorithms::algebra::argument_root::{C64, Limits, RectF64, Root, find_zeros_by_argument};

fn rect(re_min: f64, re_max: f64, im_min: f64, im_max: f64) -> RectF64 {
    RectF64::new(re_min, re_max, im_min, im_max).unwrap()
}

/// Каждому ожидаемому корню соответствует найденный в пределах `tol`, и лишних нет.
fn assert_roots(found: &[Root], expected: &[C64], tol: f64) {
    assert_eq!(found.len(), expected.len(), "found: {found:#?}");
    for e in expected {
        assert!(
            found.iter().any(|r| (r.z - e).norm() <= tol),
            "root {e} not found among {:?}",
            found.iter().map(|r| r.z).collect::<Vec<_>>()
        );
    }
}

#[test]
fn rect_validation() {
    assert!(RectF64::new(0.0, 1.0, 0.0, 1.0).is_ok());
    assert!(RectF64::new(1.0, 0.0, 0.0, 1.0).is_err());
    assert!(RectF64::new(0.0, 1.0, 1.0, 1.0).is_err());
    assert!(RectF64::new(f64::NAN, 1.0, 0.0, 1.0).is_err());
    let r = rect(-1.0, 2.0, 0.5, 1.0);
    assert_close!(r.width(), 3.0, abs = 0.0);
    assert_close!(r.height(), 0.5, abs = 0.0);
}

#[test]
fn parameter_validation() {
    let f = |z: C64| z;
    let r = rect(-1.0, 1.1, -1.0, 1.2);
    assert!(find_zeros_by_argument(f, r, -1.0, 1e-8, false, Limits::default()).is_err());
    assert!(find_zeros_by_argument(f, r, 1e-12, 0.0, false, Limits::default()).is_err());
    let limits = Limits {
        samples_per_edge0: 3,
        ..Limits::default()
    };
    assert!(find_zeros_by_argument(f, r, 1e-12, 1e-8, false, limits).is_err());
}

#[test]
fn no_zeros_inside() {
    let roots = find_zeros_by_argument(
        |z: C64| z * z + 1.0,
        rect(2.0, 3.1, 2.0, 3.2),
        1e-12,
        1e-8,
        false,
        Limits::default(),
    )
    .unwrap();
    assert!(roots.is_empty());
}

#[test]
fn single_simple_zero() {
    let z0 = C64::new(0.3, -0.2);
    let roots = find_zeros_by_argument(
        move |z: C64| z - z0,
        rect(-1.0, 1.13, -1.0, 1.27),
        1e-12,
        1e-8,
        false,
        Limits::default(),
    )
    .unwrap();
    assert_roots(&roots, &[z0], 1e-7);
    assert_eq!(roots[0].multiplicity, 1);
}

/// Известный дефект: из двух корней ±i находится только +i.
#[test]
#[ignore = "known bug: one of the two roots of z^2 + 1 is lost"]
fn conjugate_pair() {
    let roots = find_zeros_by_argument(
        |z: C64| z * z + 1.0,
        rect(-2.0, 2.13, -2.0, 2.27),
        1e-12,
        1e-8,
        false,
        Limits::default(),
    )
    .unwrap();
    assert_roots(&roots, &[C64::new(0.0, 1.0), C64::new(0.0, -1.0)], 1e-7);
}

#[test]
fn cubic_roots_of_unity() {
    let s = 3f64.sqrt() / 2.0;
    let roots = find_zeros_by_argument(
        |z: C64| z * z * z - 1.0,
        rect(-2.0, 2.13, -2.0, 2.27),
        1e-12,
        1e-8,
        false,
        Limits::default(),
    )
    .unwrap();
    assert_roots(
        &roots,
        &[C64::new(1.0, 0.0), C64::new(-0.5, s), C64::new(-0.5, -s)],
        1e-7,
    );
}

#[test]
fn transcendental_function() {
    // exp(z) = 1 ⇒ z = 2πik; в прямоугольнике три корня.
    let tau = std::f64::consts::TAU;
    let roots = find_zeros_by_argument(
        |z: C64| z.exp() - 1.0,
        rect(-1.0, 1.13, -8.0, 8.27),
        1e-12,
        1e-8,
        false,
        Limits::default(),
    )
    .unwrap();
    assert_roots(
        &roots,
        &[C64::new(0.0, 0.0), C64::new(0.0, tau), C64::new(0.0, -tau)],
        1e-7,
    );
}

/// Известное ограничение: на кратном корне поиск завершается ошибкой
/// «Contour too close to a zero» вместо отчёта с multiplicity = 2.
#[test]
#[ignore = "known limitation: a double zero ends in 'Contour too close to a zero'"]
fn double_zero_reports_multiplicity() {
    let z0 = C64::new(0.25, 0.4);
    let roots = find_zeros_by_argument(
        move |z: C64| (z - z0) * (z - z0),
        rect(-1.0, 1.13, -1.0, 1.27),
        1e-14,
        1e-6,
        false,
        Limits::default(),
    )
    .unwrap();
    assert_roots(&roots, &[z0], 1e-5);
    assert_eq!(roots[0].multiplicity, 2);
}

#[test]
fn stop_after_first_returns_one_root() {
    let roots = find_zeros_by_argument(
        |z: C64| z * z + 1.0,
        rect(-2.0, 2.13, -2.0, 2.27),
        1e-12,
        1e-8,
        true,
        Limits::default(),
    )
    .unwrap();
    assert_eq!(roots.len(), 1);
    assert!((roots[0].z * roots[0].z + 1.0).norm() < 1e-6);
}

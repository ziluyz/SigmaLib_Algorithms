mod common;

use algorithms::algebra::argument_root::{
    C64, Limits, RectF64, Root, StopReason, find_zeros_by_argument,
};

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

/// Регрессия: из двух корней ±i находился только +i (ноль рядом с линией деления давал
/// нулевое число оборотов у всех четырёх частей).
#[test]
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

/// Регрессия: на кратном корне поиск завершался ошибкой «Contour too close to a zero».
#[test]
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

/// Нули точно на линиях деления пополам (симметричный прямоугольник): линии сдвигаются.
#[test]
fn zeros_on_dividing_lines() {
    let roots = find_zeros_by_argument(
        |z: C64| z * z + 1.0,
        rect(-2.0, 2.0, -2.0, 2.0),
        1e-12,
        1e-8,
        false,
        Limits::default(),
    )
    .unwrap();
    assert_roots(&roots, &[C64::new(0.0, 1.0), C64::new(0.0, -1.0)], 1e-7);

    // z⁴ = 1: все четыре корня на осях симметрии.
    let roots = find_zeros_by_argument(
        |z: C64| z.powi(4) - 1.0,
        rect(-2.0, 2.0, -2.0, 2.0),
        1e-12,
        1e-8,
        false,
        Limits::default(),
    )
    .unwrap();
    assert_roots(
        &roots,
        &[
            C64::new(1.0, 0.0),
            C64::new(-1.0, 0.0),
            C64::new(0.0, 1.0),
            C64::new(0.0, -1.0),
        ],
        1e-7,
    );
}

/// Двойной ноль точно в центре прямоугольника: малое |f| в центре не обрывает поиск.
#[test]
fn double_zero_at_center() {
    let roots = find_zeros_by_argument(
        |z: C64| z * z,
        rect(-1.0, 1.0, -1.0, 1.0),
        1e-12,
        1e-6,
        false,
        Limits::default(),
    )
    .unwrap();
    assert_roots(&roots, &[C64::new(0.0, 0.0)], 1e-6);
    assert_eq!(roots[0].multiplicity, 2);
}

#[test]
fn triple_zero_and_simple_zero() {
    let z0 = C64::new(-0.3, 0.45);
    let z1 = C64::new(0.6, -0.1);
    let roots = find_zeros_by_argument(
        move |z: C64| (z - z0).powi(3) * (z - z1),
        rect(-1.0, 1.13, -1.0, 1.27),
        1e-14,
        1e-6,
        false,
        Limits::default(),
    )
    .unwrap();
    assert_roots(&roots, &[z0, z1], 1e-6);
    for r in &roots {
        let expected = if (r.z - z0).norm() < 1e-3 { 3 } else { 1 };
        assert_eq!(r.multiplicity, expected, "{r:?}");
    }
}

/// Многочлен с двойным корнем в развёрнутом виде: вблизи корня значения — шум округления.
/// Поиск не падает, а сообщает прямоугольник с двумя нулями.
#[test]
fn double_zero_in_rounding_noise() {
    // (z - 1)² = z² - 2z + 1
    let roots = find_zeros_by_argument(
        |z: C64| z * z - 2.0 * z + 1.0,
        rect(-1.0, 2.13, -1.0, 1.27),
        0.0,
        1e-12,
        false,
        Limits::default(),
    )
    .unwrap();
    assert_eq!(roots.iter().map(|r| r.multiplicity).sum::<i64>(), 2);
    for r in &roots {
        assert!((r.z - 1.0).norm() < 1e-6, "{r:?}");
    }
}

/// Близкие, но различные нули разделяются.
#[test]
fn close_pair_is_separated() {
    let z0 = C64::new(0.2, 0.3);
    let z1 = z0 + C64::new(1e-5, 2e-5);
    let roots = find_zeros_by_argument(
        move |z: C64| (z - z0) * (z - z1),
        rect(-1.0, 1.13, -1.0, 1.27),
        0.0,
        1e-8,
        false,
        Limits::default(),
    )
    .unwrap();
    assert_roots(&roots, &[z0, z1], 1e-7);
}

/// Число нулей с учётом кратности сохраняется: sin(z) в полосе имеет нули kπ.
#[test]
fn many_zeros_of_sine() {
    let roots = find_zeros_by_argument(
        |z: C64| z.sin(),
        rect(-10.0, 10.13, -1.0, 1.27),
        1e-12,
        1e-8,
        false,
        Limits::default(),
    )
    .unwrap();
    let expected: Vec<C64> = (-3..=3)
        .map(|k| C64::new(k as f64 * std::f64::consts::PI, 0.0))
        .collect();
    assert_roots(&roots, &expected, 1e-7);
    assert!(roots.iter().all(|r| r.multiplicity == 1));
}

/// Ноль на границе исходного прямоугольника — понятная ошибка.
#[test]
fn zero_on_boundary_is_an_error() {
    let err = find_zeros_by_argument(
        |z: C64| z * z + 1.0,
        rect(-2.0, 2.13, 1.0, 2.27),
        1e-12,
        1e-8,
        false,
        Limits::default(),
    )
    .unwrap_err();
    assert!(err.contains("boundary"), "{err}");
}

#[test]
fn stop_reasons() {
    let z0 = C64::new(0.3, -0.2);
    let search = |abs_tol: f64, rel_tol: f64, limits: Limits| {
        find_zeros_by_argument(
            move |z: C64| z - z0,
            rect(-1.0, 1.13, -1.0, 1.27),
            abs_tol,
            rel_tol,
            false,
            limits,
        )
    };

    let roots = search(1e-3, 1e-12, Limits::default()).unwrap();
    assert_eq!(roots[0].stop_reason, StopReason::AbsTolSatisfied);
    assert!(roots[0].abs_f <= 1e-3);

    let roots = search(0.0, 1e-4, Limits::default()).unwrap();
    assert_eq!(roots[0].stop_reason, StopReason::RelTolSatisfied);
    assert!(roots[0].rect_diag <= 1e-4);
    assert_roots(&roots, &[z0], 1e-4);

    let limits = Limits {
        max_depth: 5,
        ..Limits::default()
    };
    let roots = search(0.0, 1e-12, limits).unwrap();
    assert_eq!(roots.len(), 1);
    assert_eq!(roots[0].stop_reason, StopReason::LimitsReached);
    let r = roots[0].rect;
    assert!(r.re_min <= z0.re && z0.re <= r.re_max && r.im_min <= z0.im && z0.im <= r.im_max);
}

#[test]
fn limits_are_enforced() {
    let f = |z: C64| z * z * z - 1.0;
    let r = rect(-2.0, 2.13, -2.0, 2.27);

    let limits = Limits {
        max_rects: 5,
        ..Limits::default()
    };
    let err = find_zeros_by_argument(f, r, 1e-12, 1e-8, false, limits).unwrap_err();
    assert!(err.contains("budget"), "{err}");

    // Глубина, при которой индексы сетки не помещаются в u64.
    let limits = Limits {
        max_depth: 60,
        ..Limits::default()
    };
    assert!(find_zeros_by_argument(f, r, 1e-12, 1e-8, false, limits).is_err());
}

/// Полюс считается нулём отрицательной кратности.
#[test]
fn pole_has_negative_multiplicity() {
    let p = C64::new(0.1, 0.2);
    let roots = find_zeros_by_argument(
        move |z: C64| 1.0 / (z - p),
        rect(-1.0, 1.13, -1.0, 1.27),
        0.0,
        1e-8,
        false,
        Limits::default(),
    )
    .unwrap();
    assert_roots(&roots, &[p], 1e-7);
    assert_eq!(roots[0].multiplicity, -1);
}

mod common;

use algorithms::algebra::argument_root::RectF64;
use algorithms::continuation::follow;
use algorithms::num_complex::Complex64;

/// Корень z(p) = p² + i·sin p: точное решение известно, «решатель» проверяет, что оно лежит в
/// окне, и что прогноз квадратичной экстраполяцией близок к нему.
#[test]
fn follows_a_smooth_branch() {
    let exact = |p: f64| Complex64::new(p * p, p.sin());
    let params: Vec<f64> = (0..20).map(|i| 1.0 + 0.05 * f64::from(i) * (1.0 + 0.01 * f64::from(i))).collect();
    let global = RectF64::new(-10.0, 10.0, -10.0, 10.0).unwrap();
    let first = RectF64::new(0.5, 1.5, 0.5, 1.0).unwrap();
    let mut got = Vec::new();
    follow(
        &params,
        first,
        global,
        Complex64::new(1e-6, 1e-6),
        |index, p, prediction| {
            let z = exact(p);
            let w = prediction.window;
            assert!(z.re > w.re_min && z.re < w.re_max && z.im > w.im_min && z.im < w.im_max, "{index}");
            if index >= 3 {
                assert!((prediction.guess - z).norm() < 1e-3, "{index}: {}", prediction.guess);
            }
            Ok(Some(z))
        },
        |_, p, z| {
            got.push((p, z));
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(got.len(), params.len());
}

/// Неудача в узком окне расширяет его; неудача на всей области — ошибка с номером точки,
/// найденные до неё точки обработаны.
#[test]
fn widens_and_reports_failure() {
    let params = [0.0, 1.0, 2.0, 3.0];
    let global = RectF64::new(-100.0, 100.0, -1.0, 1.0).unwrap();
    let first = RectF64::new(-1.0, 1.0, -1.0, 1.0).unwrap();
    let mut widened = 0;
    let mut processed = 0;
    let error = follow(
        &params,
        first,
        global,
        Complex64::new(0.1, 0.1),
        |index, p, prediction| {
            if prediction.widened {
                widened += 1;
            }
            let w = prediction.window;
            Ok(match index {
                // На второй точке корень далеко от первой: нужно расширение окна.
                0 => Some(Complex64::new(0.0, 0.0)),
                1 => (w.re_max >= 20.0).then_some(Complex64::new(20.0, 0.0)),
                2 => Some(Complex64::new(40.0 + p - p, 0.0)),
                _ => None,
            })
        },
        |_, _, _| {
            processed += 1;
            Ok(())
        },
    )
    .unwrap_err();
    assert!(widened > 0);
    assert_eq!(processed, 3);
    assert!(error.contains("p = 3"), "{error}");
}

#[test]
fn invalid_arguments() {
    let global = RectF64::new(0.0, 1.0, 0.0, 1.0).unwrap();
    let outside = RectF64::new(2.0, 3.0, 0.0, 1.0).unwrap();
    let solve = |_: usize, _: f64, _: &algorithms::continuation::Prediction| Ok(None);
    assert!(follow(&[0.0], outside, global, Complex64::new(0.1, 0.1), solve, |_, _, _| Ok(())).is_err());
    assert!(follow(&[0.0], global, global, Complex64::new(0.0, 0.1), solve, |_, _, _| Ok(())).is_err());
    assert!(follow(&[f64::NAN], global, global, Complex64::new(0.1, 0.1), solve, |_, _, _| Ok(())).is_err());
}

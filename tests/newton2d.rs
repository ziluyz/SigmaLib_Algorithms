mod common;

use algorithms::algebra::argument_root::RectF64;
use algorithms::algebra::{Newton2dOptions, newton_2d};
use algorithms::num_complex::Complex64;

fn options() -> Newton2dOptions {
    Newton2dOptions {
        abs_tol: 1e-12,
        step_tol: Complex64::new(0.0, 0.0),
        fd_step: Complex64::new(1e-7, 1e-7),
        max_iterations: 30,
        max_halvings: 10,
    }
}

fn batch(f: impl Fn(Complex64) -> Complex64) -> impl Fn(&[Complex64]) -> Result<Vec<Complex64>, String> {
    move |zs: &[Complex64]| Ok(zs.iter().map(|&z| f(z)).collect())
}

/// Неаналитическое отображение с обращающей ориентацию производной: f = conj(z) − (1 + 2i)
/// плюс нелинейная добавка. Корень: x = 1, y = −2 (добавка в нём нулевая).
#[test]
fn non_analytic_map_with_reversed_orientation() {
    let f = |z: Complex64| {
        let (x, y) = (z.re, z.im);
        Complex64::new(x - 1.0 + 0.3 * (x - 1.0) * (y + 2.0), -y - 2.0 + 0.2 * (x - 1.0).powi(2))
    };
    let region = RectF64::new(-1.0, 3.0, -4.0, 0.0).unwrap();
    let root = newton_2d(&batch(f), Complex64::new(0.2, -1.0), region, &options())
        .unwrap()
        .unwrap();
    assert_close!(root.z.re, 1.0, abs = 1e-10);
    assert_close!(root.z.im, -2.0, abs = 1e-10);
    assert!(root.f.norm() <= 1e-12);
    assert!(root.iterations <= 8, "{root:?}");
}

/// Корень вне области не находится, точки остаются в области.
#[test]
fn stays_inside_the_region() {
    let region = RectF64::new(0.0, 1.0, 0.0, 1.0).unwrap();
    let eval = |zs: &[Complex64]| -> Result<Vec<Complex64>, String> {
        for z in zs {
            assert!(z.re >= 0.0 && z.re <= 1.0 && z.im >= 0.0 && z.im <= 1.0, "{z}");
        }
        Ok(zs.iter().map(|&z| z - Complex64::new(5.0, 0.5)).collect())
    };
    let result = newton_2d(&eval, Complex64::new(0.5, 0.5), region, &options()).unwrap();
    assert!(result.is_none());
}

/// Остановка по длине шага: функция с ошибкой округления не достигает abs_tol.
#[test]
fn stops_by_step() {
    let f = |z: Complex64| (z - Complex64::new(0.3, 0.4)) * Complex64::new(2.0, 1.0) + 1e-9;
    let mut opts = options();
    opts.abs_tol = 0.0;
    opts.step_tol = Complex64::new(1e-12, 1e-12);
    let region = RectF64::new(0.0, 1.0, 0.0, 1.0).unwrap();
    let root = newton_2d(&batch(f), Complex64::new(0.9, 0.1), region, &opts)
        .unwrap()
        .unwrap();
    assert!(root.f.norm() < 1e-14);
}

#[test]
fn invalid_arguments_and_errors() {
    let region = RectF64::new(0.0, 1.0, 0.0, 1.0).unwrap();
    let f = batch(|z| z);
    let mut opts = options();
    assert!(newton_2d(&f, Complex64::new(2.0, 0.5), region, &opts).is_err());
    opts.fd_step = Complex64::new(0.0, 1e-7);
    assert!(newton_2d(&f, Complex64::new(0.5, 0.5), region, &opts).is_err());
    let failing = |_: &[Complex64]| -> Result<Vec<Complex64>, String> { Err("boom".into()) };
    assert_eq!(
        newton_2d(&failing, Complex64::new(0.5, 0.5), region, &options()).unwrap_err(),
        "boom"
    );
    let nan = batch(|_| Complex64::new(f64::NAN, 0.0));
    assert!(newton_2d(&nan, Complex64::new(0.5, 0.5), region, &options()).unwrap().is_none());
}

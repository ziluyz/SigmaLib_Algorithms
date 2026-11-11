mod common;

use algorithms::gsl::errors::{get_error_description, turn_off_error_handler};
use algorithms::gsl::integration::{Qawo, QawoOscillator, Romberg, qng};
use std::f64::consts::PI;

/// ∫₀¹ x·sin(ωx) dx
fn x_sin(omega: f64) -> f64 {
    omega.sin() / omega.powi(2) - omega.cos() / omega
}

/// ∫₀¹ x·cos(ωx) dx
fn x_cos(omega: f64) -> f64 {
    omega.cos() / omega.powi(2) + omega.sin() / omega - 1.0 / omega.powi(2)
}

#[test]
fn qng_polynomial_and_sine() {
    turn_off_error_handler();
    let r = qng(&mut |x: f64| x * x, 0.0, 1.0, 1e-12, 1e-12).unwrap();
    assert_close!(r.result, 1.0 / 3.0, abs = 1e-12);
    assert!(r.neval > 0);

    let r = qng(&mut |x: f64| x.sin(), 0.0, PI, 1e-12, 1e-12).unwrap();
    assert_close!(r.result, 2.0, abs = 1e-11);
    assert!(r.abserr <= 1e-10);
}

#[test]
fn qng_reports_failure_as_err() {
    turn_off_error_handler();
    // Неинтегрируемая особенность: qng обязан вернуть Err, а не abort().
    let err = qng(&mut |x: f64| 1.0 / x, 0.0, 1.0, 1e-12, 1e-12).err().unwrap();
    assert!(err.contains("gsl_integration_qng"), "{err}");
}

#[test]
fn qng_captures_environment() {
    turn_off_error_handler();
    let k = 3.0;
    let mut calls = 0usize;
    let r = qng(
        &mut |x: f64| {
            calls += 1;
            k * x
        },
        0.0,
        2.0,
        1e-12,
        1e-12,
    )
    .unwrap();
    assert_close!(r.result, 6.0, abs = 1e-12);
    assert_eq!(calls, r.neval);
}

#[test]
fn romberg_requires_initialize() {
    let romberg = Romberg::new(0.0, 1.0);
    assert!(romberg.integrate(&mut |x: f64| x, 1e-10, 1e-10).is_err());
}

#[test]
fn romberg_smooth_integrands() {
    turn_off_error_handler();
    let mut romberg = Romberg::new(0.0, PI);
    romberg.initialize(20);
    let r = romberg.integrate(&mut |x: f64| x.sin(), 1e-12, 1e-12).unwrap();
    assert_close!(r.result, 2.0, abs = 1e-10);

    let mut romberg = Romberg::new(-1.0, 2.0);
    romberg.initialize(20);
    let r = romberg.integrate(&mut |x: f64| x.exp(), 1e-12, 1e-12).unwrap();
    assert_close!(r.result, 2f64.exp() - (-1f64).exp(), abs = 1e-10);
}

#[test]
fn qawo_requires_initialize() {
    let qawo = Qawo::new(0.0, 1.0);
    assert!(qawo.integrate(&mut |x: f64| x, 1e-10, 1e-10, false).is_err());
}

#[test]
fn qawo_sine_and_cosine() {
    turn_off_error_handler();
    for &omega in &[1.0, 10.0, 137.0] {
        let mut qawo = Qawo::new(0.0, 1.0);
        qawo.initialize(omega, 1000, QawoOscillator::Sine);
        let r = qawo.integrate(&mut |x: f64| x, 1e-12, 1e-10, true).unwrap();
        assert_close!(r.result, x_sin(omega), abs = 1e-10);

        let mut qawo = Qawo::new(0.0, 1.0);
        qawo.initialize(omega, 1000, QawoOscillator::Cosine);
        let r = qawo.integrate(&mut |x: f64| x, 1e-12, 1e-10, true).unwrap();
        assert_close!(r.result, x_cos(omega), abs = 1e-10);
    }
}

#[test]
fn qawo_change_parameters_and_limits() {
    turn_off_error_handler();
    let mut qawo = Qawo::new(0.0, 1.0);
    qawo.initialize(5.0, 1000, QawoOscillator::Sine);

    // другая частота и тип осциллятора на том же отрезке
    qawo.change_parameters(20.0, 0.0, 1.0, QawoOscillator::Cosine);
    let r = qawo.integrate(&mut |x: f64| x, 1e-12, 1e-10, true).unwrap();
    assert_close!(r.result, x_cos(20.0), abs = 1e-10);

    // другой отрезок: ∫₀^π cos(2x)·1 dx = 0, ∫₀^{π/4} cos(2x) dx = 1/2
    qawo.change_parameters(2.0, 0.0, PI, QawoOscillator::Cosine);
    let r = qawo.integrate(&mut |_| 1.0, 1e-12, 1e-10, true).unwrap();
    assert_close!(r.result, 0.0, abs = 1e-10);

    qawo.change_limits(0.0, PI / 4.0);
    let r = qawo.integrate(&mut |_| 1.0, 1e-12, 1e-10, true).unwrap();
    assert_close!(r.result, 0.5, abs = 1e-10);
}

#[test]
fn error_description_is_readable() {
    assert_eq!(get_error_description(0), "success");
    assert!(!get_error_description(18).is_empty());
}

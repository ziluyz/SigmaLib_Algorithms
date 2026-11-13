mod common;

use algorithms::num_complex::Complex64;
use algorithms::ode::Rk4;

/// Интегрирует y' = y от 0 до 1 за n шагов, возвращает ошибку относительно e.
fn exp_error(n: usize) -> f64 {
    let mut stepper = Rk4::<f64>::new(1);
    let mut y = [1.0];
    let mut rhs = |_x: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = y[0];
        Ok(())
    };
    for i in 0..n {
        let (x0, x1) = (i as f64 / n as f64, (i + 1) as f64 / n as f64);
        stepper.step(&mut rhs, x0, x1, &mut y).unwrap();
    }
    (y[0] - std::f64::consts::E).abs()
}

#[test]
fn rk4_has_fourth_order() {
    let (coarse, fine) = (exp_error(10), exp_error(20));
    // Удвоение числа шагов уменьшает ошибку в 2^4 = 16 раз.
    assert_close!(coarse / fine, 16.0, rel = 0.1);
    assert!(fine < 2e-7, "error {fine:e}");
}

/// Неавтономное уравнение y' = 3x², y(0) = 0: RK4 интегрирует кубический многочлен точно.
#[test]
fn rk4_uses_the_argument() {
    let mut stepper = Rk4::<f64>::new(1);
    let mut y = [0.0];
    stepper
        .step(
            &mut |x, _y: &[f64], dy: &mut [f64]| {
                dy[0] = 3.0 * x * x;
                Ok(())
            },
            0.0,
            2.0,
            &mut y,
        )
        .unwrap();
    assert_close!(y[0], 8.0, abs = 1e-13);
}

/// Комплексное состояние: y' = i·y, решение exp(i·x); шаг назад возвращает в начало.
#[test]
fn rk4_complex_state_and_backward_step() {
    let mut stepper = Rk4::<Complex64>::new(1);
    let mut y = [Complex64::ONE];
    let mut rhs = |_x: f64, y: &[Complex64], dy: &mut [Complex64]| {
        dy[0] = Complex64::I * y[0];
        Ok(())
    };
    let n = 100;
    for i in 0..n {
        stepper
            .step(&mut rhs, i as f64 / n as f64, (i + 1) as f64 / n as f64, &mut y)
            .unwrap();
    }
    assert_close!(y[0].re, 1f64.cos(), abs = 1e-9);
    assert_close!(y[0].im, 1f64.sin(), abs = 1e-9);
    for i in (0..n).rev() {
        stepper
            .step(&mut rhs, (i + 1) as f64 / n as f64, i as f64 / n as f64, &mut y)
            .unwrap();
    }
    assert_close!(y[0].re, 1.0, abs = 1e-9);
    assert_close!(y[0].im, 0.0, abs = 1e-9);
}

#[test]
fn rk4_errors() {
    let mut stepper = Rk4::<f64>::new(2);
    let mut short = [0.0];
    assert!(
        stepper
            .step(&mut |_, _: &[f64], _: &mut [f64]| Ok(()), 0.0, 1.0, &mut short)
            .is_err()
    );
    let mut y = [1.0, 2.0];
    let error = stepper
        .step(
            &mut |_, _: &[f64], _: &mut [f64]| Err("boom".to_string()),
            0.0,
            1.0,
            &mut y,
        )
        .unwrap_err();
    assert_eq!(error, "boom");
    assert_eq!(y, [1.0, 2.0]);
}

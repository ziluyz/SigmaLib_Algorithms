mod common;

use algorithms::algebra::find_root_sequence;
use algorithms::gsl::errors::turn_off_error_handler;
use algorithms::gsl::sfunctions::bessel::{djn, jn};

// Эталоны: Abramowitz & Stegun, табл. 9.1 и 9.5; DLMF §10.21.
#[test]
fn jn_reference_values() {
    assert_close!(jn(0, 0.0).unwrap(), 1.0, abs = 1e-15);
    assert_close!(jn(1, 0.0).unwrap(), 0.0, abs = 1e-15);
    assert_close!(jn(0, 1.0).unwrap(), 0.765_197_686_557_966_6, abs = 1e-14);
    assert_close!(jn(1, 1.0).unwrap(), 0.440_050_585_744_933_5, abs = 1e-14);
    assert_close!(jn(2, 5.0).unwrap(), 0.046_565_116_277_752_2, abs = 1e-14);
}

#[test]
fn jn_negative_order_symmetry() {
    for n in 0..6 {
        for &x in &[0.3, 1.7, 4.2, 11.0] {
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            assert_close!(jn(-n, x).unwrap(), sign * jn(n, x).unwrap(), abs = 1e-14);
        }
    }
}

#[test]
fn jn_recurrence() {
    // J_{n-1}(x) + J_{n+1}(x) = (2n/x) J_n(x)
    for n in 1..8 {
        for &x in &[0.5, 2.0, 7.5, 20.0] {
            let lhs = jn(n - 1, x).unwrap() + jn(n + 1, x).unwrap();
            let rhs = 2.0 * n as f64 / x * jn(n, x).unwrap();
            assert_close!(lhs, rhs, abs = 1e-13);
        }
    }
}

#[test]
fn djn_matches_identities() {
    for &x in &[0.4, 1.3, 6.1] {
        // J0' = -J1
        assert_close!(djn(0, x).unwrap(), -jn(1, x).unwrap(), abs = 1e-14);
        // центральная разность
        let h = 1e-5;
        for n in 0..5 {
            let numeric = (jn(n, x + h).unwrap() - jn(n, x - h).unwrap()) / (2.0 * h);
            assert_close!(djn(n, x).unwrap(), numeric, abs = 1e-8);
        }
    }
}

#[test]
fn derivative_zeros_match_tables() {
    // Нули J_m'(x) — так их ищут проекты (см. Sigma_GyrStCur).
    let roots = find_root_sequence(&|x| djn(1, x), 1.1, 3.0, 3).unwrap();
    assert_close!(roots[0], 1.841_183_781_340_659, abs = 1e-6);
    assert_close!(roots[1], 5.331_442_773_525_033, abs = 1e-6);
    assert_close!(roots[2], 8.536_316_366_346_286, abs = 1e-6);

    let roots = find_root_sequence(&|x| djn(0, x), 0.1, 3.0, 2).unwrap();
    assert_close!(roots[0], 3.831_705_970_207_512, abs = 1e-6);
    assert_close!(roots[1], 7.015_586_669_815_619, abs = 1e-6);
}

/// Исчезающе малое значение — не ошибка: GSL сообщает underflow, обёртка возвращает ноль.
#[test]
fn jn_underflow_is_zero() {
    turn_off_error_handler();
    assert_close!(jn(1000, 1e-3).unwrap(), 0.0, abs = 1e-300);
    assert_close!(djn(1000, 1e-3).unwrap(), 0.0, abs = 1e-300);
}

/// Нечисловой аргумент — `Err`, а не молча NaN.
#[test]
fn jn_reports_gsl_error() {
    turn_off_error_handler();
    match jn(1, f64::NAN) {
        Err(err) => assert!(err.contains("gsl_sf_bessel_Jn"), "{err}"),
        Ok(v) => panic!("expected an error, got {v}"),
    }
}

#[test]
fn series_matches_gsl() {
    use algorithms::gsl::sfunctions::bessel::{SERIES_LIMIT, djn, jn, jn_djn};
    algorithms::gsl::errors::turn_off_error_handler();
    for n in -6..=6 {
        for i in -40..=40 {
            let x = SERIES_LIMIT * f64::from(i) / 40.0;
            let (j, dj) = jn_djn(n, x).unwrap();
            assert_close!(j, jn(n, x).unwrap(), abs = 1e-16, rel = 1e-14);
            // djn у GSL — разность J_{n−1} и J_{n+1}: абсолютная ошибка около 1e-16.
            assert_close!(dj, djn(n, x).unwrap(), abs = 1e-15, rel = 1e-14);
        }
    }
    // За пределом ряда — значения GSL.
    let (j, dj) = jn_djn(3, 7.5).unwrap();
    assert_eq!(j, jn(3, 7.5).unwrap());
    assert_eq!(dj, djn(3, 7.5).unwrap());
    // Табличные значения (Abramowitz & Stegun, табл. 9.1): J_0(1), J_1(1); J_1'(0) = 1/2.
    assert_close!(jn_djn(0, 1.0).unwrap().0, 0.765_197_686_557_966_6, rel = 1e-15);
    assert_close!(jn_djn(0, 1.0).unwrap().1, -0.440_050_585_744_933_5, rel = 2e-15);
    assert_eq!(jn_djn(1, 0.0).unwrap(), (0.0, 0.5));
    assert!(jn_djn(1, f64::NAN).is_err());
}

#[test]
fn bessel_over_power() {
    use algorithms::gsl::sfunctions::bessel::{jn, jn_over_power};
    algorithms::gsl::errors::turn_off_error_handler();
    for m in 0..=5u32 {
        // В нуле — 1/(2^m m!).
        let factorial: f64 = (1..=m).map(f64::from).product();
        assert_close!(jn_over_power(m, 0.0).unwrap(), 1.0 / (2f64.powi(m as i32) * factorial), rel = 1e-15);
        for x in [-3.0, -1.0, 1e-3, 0.5, 1.9, 2.1, 4.0] {
            let expected = jn(m as i32, x).unwrap() / x.powi(m as i32);
            assert_close!(jn_over_power(m, x).unwrap(), expected, abs = 1e-16, rel = 1e-13);
        }
    }
    assert!(jn_over_power(1, f64::INFINITY).is_err());
}

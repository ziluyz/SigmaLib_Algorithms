mod common;

use algorithms::linalg::{CMatrix, RMatrix, expm, solve};
use algorithms::num_complex::Complex64;

fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

#[test]
fn solve_complex_system() {
    let a = CMatrix::from_row_slice(2, 2, &[c(1.0, 1.0), c(2.0, 0.0), c(0.0, -1.0), c(3.0, 1.0)]);
    let x = CMatrix::from_row_slice(2, 2, &[c(1.0, 0.0), c(0.0, 2.0), c(-1.0, 1.0), c(0.5, 0.0)]);
    let b = &a * &x;
    let found = solve(&a, &b).unwrap();
    for (found, expected) in found.iter().zip(x.iter()) {
        assert_close!(found.re, expected.re, abs = 1e-13);
        assert_close!(found.im, expected.im, abs = 1e-13);
    }
}

#[test]
fn solve_errors() {
    let singular = CMatrix::from_row_slice(2, 2, &[c(1.0, 0.0), c(2.0, 0.0), c(2.0, 0.0), c(4.0, 0.0)]);
    let b = CMatrix::from_row_slice(2, 1, &[c(1.0, 0.0), c(1.0, 0.0)]);
    assert!(solve(&singular, &b).is_err());
    let zero = CMatrix::zeros(2, 2);
    assert!(solve(&zero, &b).is_err());
    let wide = CMatrix::zeros(2, 3);
    assert!(solve(&wide, &b).is_err());
}

/// exp([[0, -t], [t, 0]]) — поворот на угол t.
#[test]
fn expm_rotation() {
    let t = 2.5;
    let e = expm(&RMatrix::from_row_slice(2, 2, &[0.0, -t, t, 0.0])).unwrap();
    assert_close!(e[(0, 0)], t.cos(), abs = 1e-13);
    assert_close!(e[(0, 1)], -t.sin(), abs = 1e-13);
    assert_close!(e[(1, 0)], t.sin(), abs = 1e-13);
    assert_close!(e[(1, 1)], t.cos(), abs = 1e-13);
}

/// Нильпотентная матрица: ряд обрывается, exp(N) = 1 + N + N²/2.
#[test]
fn expm_nilpotent_and_large_norm() {
    let n = RMatrix::from_row_slice(3, 3, &[0.0, 40.0, 7.0, 0.0, 0.0, 30.0, 0.0, 0.0, 0.0]);
    let e = expm(&n).unwrap();
    assert_close!(e[(0, 1)], 40.0, rel = 1e-12);
    assert_close!(e[(1, 2)], 30.0, rel = 1e-12);
    assert_close!(e[(0, 2)], 7.0 + 600.0, rel = 1e-12);
    assert_close!(e[(2, 2)], 1.0, rel = 1e-12);
    assert!(expm(&RMatrix::zeros(2, 3)).is_err());
    assert_eq!(expm(&RMatrix::zeros(0, 0)).unwrap().nrows(), 0);
}

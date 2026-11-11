mod common;

use algorithms::algebra::find_root_sequence;
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

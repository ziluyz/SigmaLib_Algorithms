//! Общие помощники интеграционных тестов.

/// Сравнение чисел с явным допуском: `|a - b| <= abs + rel * max(|a|, |b|)`.
#[macro_export]
macro_rules! assert_close {
    ($a:expr, $b:expr, abs = $abs:expr) => {
        $crate::assert_close!($a, $b, abs = $abs, rel = 0.0)
    };
    ($a:expr, $b:expr, rel = $rel:expr) => {
        $crate::assert_close!($a, $b, abs = 0.0, rel = $rel)
    };
    ($a:expr, $b:expr, abs = $abs:expr, rel = $rel:expr) => {{
        let (a, b): (f64, f64) = ($a, $b);
        let tol = $abs + $rel * a.abs().max(b.abs());
        assert!(
            (a - b).abs() <= tol,
            "assert_close failed: {} = {:e}, {} = {:e}, |diff| = {:e} > tol = {:e}",
            stringify!($a),
            a,
            stringify!($b),
            b,
            (a - b).abs(),
            tol
        );
    }};
}

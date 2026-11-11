mod common;

use algorithms::algebra::{find_root, find_root_sequence};
use std::f64::consts::PI;

#[test]
fn find_root_sqrt2() {
    let x = find_root(&|x| Ok(x * x - 2.0), 0.0, 2.0, 1e-12).unwrap();
    assert_close!(x, 2f64.sqrt(), abs = 1e-12);
}

#[test]
fn find_root_decreasing_function() {
    let x = find_root(&|x| Ok(1.0 - x), -3.0, 5.5, 1e-10).unwrap();
    assert_close!(x, 1.0, abs = 1e-10);
}

/// Регрессия: середина отрезка точно попадает в корень (раньше возвращался x_max = 5).
#[test]
fn find_root_midpoint_hits_root_exactly() {
    let x = find_root(&|x| Ok(1.0 - x), -3.0, 5.0, 1e-10).unwrap();
    assert_close!(x, 1.0, abs = 1e-10);
}

/// Регрессия: корень точно на левой границе (раньше возвращался x_max).
#[test]
fn find_root_at_left_border() {
    let x = find_root(&|x| Ok(x), 0.0, 4.0, 1e-10).unwrap();
    assert_close!(x, 0.0, abs = 1e-10);
}

#[test]
fn find_root_at_right_border() {
    let x = find_root(&|x| Ok(x - 4.0), 0.0, 4.0, 1e-10).unwrap();
    assert_close!(x, 4.0, abs = 1e-10);
}

#[test]
fn find_root_rejects_same_sign() {
    let err = find_root(&|x| Ok(x * x + 1.0), -1.0, 1.0, 1e-9).unwrap_err();
    assert!(err.contains("same sign"), "{err}");
}

#[test]
fn find_root_rejects_inverted_interval() {
    assert!(find_root(&|x| Ok(x), 1.0, -1.0, 1e-9).is_err());
    assert!(find_root(&|x| Ok(x), 1.0, 1.0, 1e-9).is_err());
}

#[test]
fn find_root_propagates_function_error() {
    let err = find_root(&|_| Err("boom".to_string()), 0.0, 1.0, 1e-9).unwrap_err();
    assert_eq!(err, "boom");
}

#[test]
fn find_root_sequence_zeros_of_sine() {
    let roots = find_root_sequence(&|x: f64| Ok(x.sin()), 0.1, 3.0, 5).unwrap();
    assert_eq!(roots.len(), 5);
    for (k, root) in roots.iter().enumerate() {
        assert_close!(*root, PI * (k + 1) as f64, abs = 1e-6);
    }
}

#[test]
fn find_root_sequence_zero_count() {
    let roots = find_root_sequence(&|x: f64| Ok(x.sin()), 0.1, 3.0, 0).unwrap();
    assert!(roots.is_empty());
}

/// Известный дефект (ROADMAP): если корней меньше, чем запрошено, поиск не завершается.
#[test]
#[ignore = "known bug: find_root_sequence never returns when fewer than n roots exist"]
fn find_root_sequence_terminates_when_roots_run_out() {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        // У x - 1 единственный корень.
        let _ = tx.send(find_root_sequence(&|x| Ok(x - 1.0), 0.0, 1.0, 2));
    });
    let result = rx
        .recv_timeout(std::time::Duration::from_secs(2))
        .expect("find_root_sequence did not return within 2 s");
    assert!(result.is_err());
}

/// Известный дефект (ROADMAP): ошибка функции неотличима от «нет корня на отрезке».
#[test]
#[ignore = "known bug: find_root_sequence swallows function errors and loops forever"]
fn find_root_sequence_propagates_function_error() {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(find_root_sequence(&|_| Err("boom".to_string()), 0.0, 1.0, 1));
    });
    let result = rx
        .recv_timeout(std::time::Duration::from_secs(2))
        .expect("find_root_sequence did not return within 2 s");
    assert_eq!(result.unwrap_err(), "boom");
}

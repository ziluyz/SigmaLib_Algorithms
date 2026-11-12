mod common;

use algorithms::algebra::{
    MAX_ROOT_GAP, find_root, find_root_sequence, find_root_sequence_bounded,
};
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

/// Регрессия: если корней меньше, чем запрошено, поиск не завершался.
#[test]
fn find_root_sequence_terminates_when_roots_run_out() {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        // У x - 1 единственный корень.
        let _ = tx.send(find_root_sequence(&|x| Ok(x - 1.0), 0.0, 1.0, 2));
    });
    let result = rx
        .recv_timeout(std::time::Duration::from_secs(2))
        .expect("find_root_sequence did not return within 2 s");
    let err = result.unwrap_err();
    assert!(err.contains("found 1 roots of 2"), "{err}");
}

/// Регрессия: ошибка функции была неотличима от «нет корня на отрезке», поиск не завершался.
#[test]
fn find_root_sequence_propagates_function_error() {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(find_root_sequence(
            &|_| Err("boom".to_string()),
            0.0,
            1.0,
            1,
        ));
    });
    let result = rx
        .recv_timeout(std::time::Duration::from_secs(2))
        .expect("find_root_sequence did not return within 2 s");
    assert_eq!(result.unwrap_err(), "boom");
}

/// Ошибка функции вдали от стартовой точки тоже доходит до вызывающего.
#[test]
fn find_root_sequence_propagates_late_function_error() {
    let f = |x: f64| {
        if x > 7.0 {
            Err(format!("out of domain: {x}"))
        } else {
            Ok(x.sin())
        }
    };
    let err = find_root_sequence(&f, 0.1, 3.0, 5).unwrap_err();
    assert!(err.starts_with("out of domain"), "{err}");
    // Двух корней до границы области хватает.
    assert_eq!(find_root_sequence(&f, 0.1, 3.0, 2).unwrap().len(), 2);
}

/// Корень дальше предела поиска — ошибка, а не бесконечный цикл.
#[test]
fn find_root_sequence_search_limit() {
    let far = 2.0 * MAX_ROOT_GAP;
    let err = find_root_sequence(&|x| Ok(x - far), 0.0, 1.0, 1).unwrap_err();
    assert!(err.contains("found 0 roots of 1"), "{err}");

    let near = 0.5 * MAX_ROOT_GAP;
    let roots = find_root_sequence(&|x| Ok(x - near), 0.0, 1.0, 1).unwrap();
    assert_close!(roots[0], near, abs = 1e-6);
}

#[test]
fn find_root_sequence_rejects_bad_arguments() {
    let f = |x: f64| Ok(x.sin());
    assert!(find_root_sequence(&f, 0.1, 0.0, 1).is_err());
    assert!(find_root_sequence(&f, 0.1, -1.0, 1).is_err());
    assert!(find_root_sequence(&f, 0.1, f64::NAN, 1).is_err());
    assert!(find_root_sequence(&f, f64::INFINITY, 1.0, 1).is_err());
}

/// Корень точно в стартовой точке входит в последовательность.
#[test]
fn find_root_sequence_root_at_start() {
    let roots = find_root_sequence(&|x: f64| Ok(x * (x - 2.0)), 0.0, 1.0, 2).unwrap();
    assert_close!(roots[0], 0.0, abs = 0.0);
    assert_close!(roots[1], 2.0, abs = 1e-6);
}

#[test]
fn find_root_sequence_bounded_stops_at_border() {
    let f = |x: f64| Ok(x.sin());
    // На [0.1, 10] три нуля синуса: π, 2π, 3π.
    let roots = find_root_sequence_bounded(&f, 0.1, 10.0, 3.0, 100).unwrap();
    assert_eq!(roots.len(), 3);
    for (k, root) in roots.iter().enumerate() {
        assert_close!(*root, PI * (k + 1) as f64, abs = 1e-6);
    }
    // Не больше запрошенного числа.
    assert_eq!(
        find_root_sequence_bounded(&f, 0.1, 10.0, 3.0, 2)
            .unwrap()
            .len(),
        2
    );
    // Корней нет — пустой результат, не ошибка.
    assert!(
        find_root_sequence_bounded(&f, 0.1, 3.0, 3.0, 5)
            .unwrap()
            .is_empty()
    );
    // Функция не вычисляется правее границы.
    let guarded = |x: f64| {
        if x > 4.0 {
            Err("beyond".to_string())
        } else {
            Ok(x.sin())
        }
    };
    let roots = find_root_sequence_bounded(&guarded, 0.1, 4.0, 3.0, 5).unwrap();
    assert_eq!(roots.len(), 1);
    assert_close!(roots[0], PI, abs = 1e-6);
}

#[test]
fn find_root_sequence_bounded_rejects_bad_border() {
    let f = |x: f64| Ok(x.sin());
    assert!(find_root_sequence_bounded(&f, 1.0, 0.0, 1.0, 1).is_err());
    assert!(find_root_sequence_bounded(&f, 1.0, f64::INFINITY, 1.0, 1).is_err());
    assert!(find_root_sequence_bounded(&f, 1.0, f64::NAN, 1.0, 1).is_err());
}

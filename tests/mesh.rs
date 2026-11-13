mod common;

use algorithms::mesh::{Interval, refine};

#[test]
fn segments_are_divided_evenly() {
    let mesh = refine(&[0.0, 1.0, 1.5], &[&[0.0, 2.0, 2.0]], 0.25).unwrap();
    assert_eq!(mesh.x, vec![0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]);
    assert_eq!(mesh.y[0], vec![0.0, 0.5, 1.0, 1.5, 2.0, 2.0, 2.0]);
    assert!(mesh.intervals.iter().all(|&i| i == Interval::Smooth));
    assert_eq!(mesh.intervals.len(), mesh.len() - 1);
    assert_close!(mesh.slope(0, 0).unwrap(), 2.0, abs = 1e-12);
    assert_close!(mesh.slope(0, 5).unwrap(), 0.0, abs = 1e-12);
}

#[test]
fn step_never_exceeds_the_limit() {
    let mesh = refine(&[0.0, 1.0], &[&[1.0, 3.0]], 0.3).unwrap();
    assert_eq!(mesh.len(), 5); // 4 шага по 0.25
    for pair in mesh.x.windows(2) {
        assert!(pair[1] - pair[0] <= 0.3);
    }
    assert_eq!(*mesh.x.last().unwrap(), 1.0);
    assert_eq!(*mesh.y[0].last().unwrap(), 3.0);
}

/// Профиль из input.ixml проекта SParameter: скачок радиуса в z = 0.
#[test]
fn jump_keeps_both_values() {
    let mesh = refine(
        &[-1.0, 0.0, 0.0, 1.0],
        &[&[0.25, 0.25, 0.3054, 0.3054]],
        0.5,
    )
    .unwrap();
    assert_eq!(mesh.x, vec![-1.0, -0.5, 0.0, 0.0, 0.5, 1.0]);
    assert_eq!(mesh.y[0], vec![0.25, 0.25, 0.25, 0.3054, 0.3054, 0.3054]);
    assert_eq!(mesh.intervals[2], Interval::Jump);
    assert_eq!(
        mesh.intervals.iter().filter(|&&i| i == Interval::Jump).count(),
        1
    );
    assert_eq!(mesh.slope(0, 2), None);
}

#[test]
fn several_channels_share_the_grid() {
    // Скачок только во втором канале; повтор точки без скачка отбрасывается.
    let mesh = refine(
        &[0.0, 1.0, 1.0, 2.0, 2.0, 3.0],
        &[&[0.0, 1.0, 1.0, 2.0, 2.0, 3.0], &[5.0, 5.0, 7.0, 7.0, 7.0, 7.0]],
        1.0,
    )
    .unwrap();
    assert_eq!(mesh.x, vec![0.0, 1.0, 1.0, 2.0, 3.0]);
    assert_eq!(mesh.y[0], vec![0.0, 1.0, 1.0, 2.0, 3.0]);
    assert_eq!(mesh.y[1], vec![5.0, 5.0, 7.0, 7.0, 7.0]);
}

#[test]
fn invalid_input_is_error() {
    assert!(refine(&[0.0], &[&[1.0]], 0.1).is_err());
    assert!(refine(&[0.0, 1.0], &[&[1.0]], 0.1).is_err());
    assert!(refine(&[0.0, 1.0], &[&[1.0, 2.0]], 0.0).is_err());
    assert!(refine(&[0.0, 1.0], &[&[1.0, 2.0]], f64::NAN).is_err());
    assert!(refine(&[0.0, 1.0, 0.5], &[&[1.0, 2.0, 3.0]], 0.1).is_err());
    assert!(refine(&[0.0, f64::NAN], &[&[1.0, 2.0]], 0.1).is_err());
    assert!(refine(&[0.0, 1.0], &[&[1.0, f64::INFINITY]], 0.1).is_err());
    assert!(refine(&[1.0, 1.0], &[&[1.0, 2.0]], 0.1).is_err());
    assert!(refine(&[0.0, 1.0], &[&[1.0, 2.0]], 1e-12).is_err());
}

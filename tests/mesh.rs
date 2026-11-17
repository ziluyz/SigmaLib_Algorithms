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

#[test]
fn limits_inside_and_at_jump() {
    use algorithms::mesh::limits;
    let xs = [0.0, 1.0, 1.0, 2.0];
    let ys = [0.0, 1.0, 3.0, 5.0];
    assert_eq!(limits(&xs, &ys, 0.5).unwrap(), (0.5, 0.5));
    assert_eq!(limits(&xs, &ys, 1.0).unwrap(), (1.0, 3.0));
    assert_eq!(limits(&xs, &ys, 0.0).unwrap(), (0.0, 0.0));
    assert_eq!(limits(&xs, &ys, 2.0).unwrap(), (5.0, 5.0));
    assert_close!(limits(&xs, &ys, 1.5).unwrap().0, 4.0, abs = 1e-15);
    assert!(limits(&xs, &ys, 2.5).is_err());
    assert!(limits(&xs, &ys[..2], 0.5).is_err());
}

/// Профили из input.ixml проекта ssmgnfs: геометрия и магнитное поле на разных сетках.
#[test]
fn union_of_two_profiles() {
    use algorithms::mesh::refine_union;
    let zg = [-0.1, 0.0, 1.2, 1.4, 2.0];
    let rg = [0.145, 0.15, 0.15, 0.16, 0.16];
    let zb = [-10.0, 1.5, 1.7, 40.0];
    let wb = [1.0, 1.0, 0.2, 0.2];
    let mesh = refine_union(&[(&zg, &rg), (&zb, &wb)], -0.1, 2.0, 0.01).unwrap();
    assert_eq!(mesh.len(), 211);
    assert!(mesh.intervals.iter().all(|&i| i == Interval::Smooth));
    for (i, &z) in mesh.x.iter().enumerate() {
        assert_close!(z, -0.1 + 0.01 * i as f64, abs = 1e-12);
    }
    // Узлы обоих профилей сохранены точно.
    for knot in [0.0, 1.2, 1.4, 1.5, 1.7] {
        assert!(mesh.x.contains(&knot), "{knot}");
    }
    let at = |z: f64| mesh.x.iter().position(|&x| x == z).unwrap();
    assert_eq!(mesh.y[0][at(1.4)], 0.16);
    assert_eq!(mesh.y[1][at(1.5)], 1.0);
    assert_close!(mesh.y[1][at(1.7) - 10], 0.6, abs = 1e-12);
    assert_close!(mesh.y[0][at(1.2) + 10], 0.155, abs = 1e-12);
}

#[test]
fn union_keeps_jumps_and_checks_coverage() {
    use algorithms::mesh::refine_union;
    let mesh = refine_union(
        &[(&[0.0, 1.0, 1.0, 2.0], &[1.0, 1.0, 2.0, 2.0]), (&[0.0, 2.0], &[0.0, 4.0])],
        0.0,
        2.0,
        0.5,
    )
    .unwrap();
    assert_eq!(mesh.x, vec![0.0, 0.5, 1.0, 1.0, 1.5, 2.0]);
    assert_eq!(mesh.y[0], vec![1.0, 1.0, 1.0, 2.0, 2.0, 2.0]);
    assert_eq!(mesh.y[1], vec![0.0, 1.0, 2.0, 2.0, 3.0, 4.0]);
    assert_eq!(mesh.intervals[2], Interval::Jump);
    assert!(refine_union(&[(&[0.5, 2.0], &[0.0, 1.0])], 0.0, 2.0, 0.5).is_err());
    assert!(refine_union(&[], 0.0, 2.0, 0.5).is_err());
    assert!(refine_union(&[(&[0.0, 2.0], &[0.0, 1.0])], 1.0, 1.0, 0.5).is_err());
}

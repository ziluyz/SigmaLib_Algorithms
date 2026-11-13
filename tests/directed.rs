mod common;

use algorithms::bvp::directed::{Chain, Policy, Transform, combine, solve_end_condition};
use algorithms::num_complex::Complex64;

/// Детерминированные «случайные» комплексные числа.
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }
    fn complex(&mut self, n: usize) -> Vec<Complex64> {
        (0..n).map(|_| Complex64::new(self.next(), self.next())).collect()
    }
}

/// `w` (r × rows) · `basis` (rows × cols), построчно.
fn product(w: &[Complex64], basis: &[Complex64], r: usize, rows: usize, cols: usize) -> Vec<Complex64> {
    let mut out = vec![Complex64::ZERO; r * cols];
    for s in 0..r {
        for j in 0..cols {
            out[s * cols + j] = (0..rows).map(|i| w[s * rows + i] * basis[i * cols + j]).sum();
        }
    }
    out
}

fn assert_vectors_close(a: &[Complex64], b: &[Complex64], tol: f64) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_close!(a.re, b.re, abs = tol);
        assert_close!(a.im, b.im, abs = tol);
    }
}

const R: usize = 3;
const ROWS: usize = 6;

/// После переразложения в m-м векторе одна растущая компонента (заданной величины), в
/// частном решении — ни одной; само решение не меняется. `cols` = R (однородная) и R + 1.
fn check_transform(cols: usize) {
    let mut rng = Lcg(cols as u64);
    let basis = rng.complex(ROWS * cols);
    let w = rng.complex(R * ROWS);
    let target = rng.complex(R);
    let growing = product(&w, &basis, R, ROWS, cols);
    let transform = Transform::new(&growing, R, cols, &target).unwrap();

    let mut new_basis = basis.clone();
    transform.apply(&mut new_basis, cols).unwrap();
    let new_growing = product(&w, &new_basis, R, ROWS, cols);
    for s in 0..R {
        for j in 0..cols {
            let expected = if s == j { target[s] } else { Complex64::ZERO };
            assert_vectors_close(&[new_growing[s * cols + j]], &[expected], 1e-11);
        }
    }

    let c_after = rng.complex(R);
    let c_before = transform.coefficients_before(&c_after).unwrap();
    let (mut after, mut before) = (vec![Complex64::ZERO; ROWS], vec![Complex64::ZERO; ROWS]);
    combine(&new_basis, cols, &c_after, &mut after).unwrap();
    combine(&basis, cols, &c_before, &mut before).unwrap();
    assert_vectors_close(&after, &before, 1e-11);
}

#[test]
fn transform_homogeneous() {
    check_transform(R);
}

/// Режим с частным решением (число столбцов R + 1) — тот, в котором ошибался оригинал
/// SParameter на C++ (матрица O читалась с шагом строки R вместо R + 1).
#[test]
fn transform_with_particular_solution() {
    check_transform(R + 1);
}

/// Цепочка: базис «прогоняется» (умножается слева на матрицу шага) и дважды переразлагается;
/// решение, восстановленное на каждом участке, — образ одного и того же начального решения.
#[test]
fn chain_restores_the_same_solution_on_every_segment() {
    let cols = R + 1;
    let mut rng = Lcg(42);
    let step: Vec<Complex64> = rng.complex(ROWS * ROWS);
    let propagate = |basis: &[Complex64]| product(&step, basis, ROWS, ROWS, cols);
    let w = rng.complex(R * ROWS);
    let target = vec![Complex64::ONE; R];

    let mut chain = Chain::new();
    let mut stored = vec![rng.complex(ROWS * cols)];
    for node in 1..=4 {
        let mut basis = propagate(&stored[node - 1]);
        // Частное решение неоднородной задачи: к последнему столбцу ничего не добавляем —
        // «шаг» линеен по всему базису, что соответствует f = 0 с b ≠ 0.
        if node == 2 || node == 3 {
            let growing = product(&w, &basis, R, ROWS, cols);
            let transform = Transform::new(&growing, R, cols, &target).unwrap();
            transform.apply(&mut basis, cols).unwrap();
            chain.push(node, transform).unwrap();
        }
        stored.push(basis);
    }
    assert_eq!(chain.len(), 2);

    let c_end = rng.complex(R);
    let segments = chain.segments(&c_end).unwrap();
    assert_eq!(
        segments.iter().map(|s| s.0).collect::<Vec<_>>(),
        vec![0, 2, 3]
    );
    let solution_at = |node: usize| {
        let c = &segments.iter().rev().find(|s| s.0 <= node).unwrap().1;
        let mut y = vec![Complex64::ZERO; ROWS];
        combine(&stored[node], cols, c, &mut y).unwrap();
        y
    };
    // y(node) = step · y(node − 1). Частное решение входит с коэффициентом 1, поэтому шаг
    // применяем к матрице из одного столбца.
    for node in 1..=4 {
        let previous = solution_at(node - 1);
        let expected = product(&step, &previous, ROWS, ROWS, 1);
        assert_vectors_close(&solution_at(node), &expected, 1e-9);
    }
    assert!(chain.push(3, Transform::new(&[], 0, 0, &[]).unwrap()).is_err());
}

#[test]
fn end_condition_is_satisfied() {
    for cols in [R, R + 1] {
        let mut rng = Lcg(7 + cols as u64);
        let basis = rng.complex(ROWS * cols);
        let d_operator = rng.complex(R * ROWS);
        let d = rng.complex(R);
        let condition = product(&d_operator, &basis, R, ROWS, cols);
        let c = solve_end_condition(&condition, R, cols, &d).unwrap();
        let mut y = vec![Complex64::ZERO; ROWS];
        combine(&basis, cols, &c, &mut y).unwrap();
        assert_vectors_close(&product(&d_operator, &y, R, ROWS, 1), &d, 1e-11);
    }
}

#[test]
fn policy() {
    let small = [Complex64::new(1.0, -2.0)];
    let large = [Complex64::new(1.0, -2e5)];
    assert!(!Policy::Never.due(3, &large));
    assert!(!Policy::Threshold(1e4).due(3, &small));
    assert!(Policy::Threshold(1e4).due(3, &large));
    assert!(Policy::Nodes(vec![3, 9]).due(3, &small));
    assert!(!Policy::Nodes(vec![3, 9]).due(4, &large));
}

#[test]
fn errors() {
    let z = Complex64::ZERO;
    // Вырожденная матрица растущих компонент.
    assert!(Transform::new(&[z; 4], 2, 2, &[Complex64::ONE; 2]).is_err());
    // Неверные размеры.
    assert!(Transform::new(&[z; 4], 2, 4, &[z; 2]).is_err());
    assert!(Transform::new(&[z; 5], 2, 2, &[z; 2]).is_err());
    let identity = [Complex64::ONE, z, z, Complex64::ONE];
    let transform = Transform::new(&identity, 2, 2, &[Complex64::ONE; 2]).unwrap();
    assert!(transform.apply(&mut [z; 5], 2).is_err());
    assert!(transform.apply(&mut [z; 8], 4).is_err());
    assert!(transform.coefficients_before(&[z; 3]).is_err());
    assert!(combine(&[z; 4], 2, &[z; 2], &mut [z; 3]).is_err());
    assert!(solve_end_condition(&[z; 4], 2, 2, &[z; 2]).is_err());
}

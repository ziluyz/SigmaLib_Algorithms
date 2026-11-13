//! Прогонка с направленной ортогонализацией (А. А. Быков, А. С. Ильинский, ЖВМиМФ, 1979,
//! т. 19, № 3, с. 631) — алгебраическая часть метода, не зависящая от уравнений.
//!
//! Краевая задача `dy/dt = A(t) y + f(t)`, `B y(0) = b`, `D y(T) = d` для вектора `y` длины
//! `l` решается прогонкой слева направо **базиса**: `r` линейно независимых решений
//! однородной системы, удовлетворяющих `B y = 0`, и (если задача неоднородна) частного
//! решения `y₀` с `B y₀ = b`. Решение — `y = y₀ + Σ cₘ yₘ`, коэффициенты `c` находятся из
//! условия на правом конце ([`solve_end_condition`]).
//!
//! Экспоненциально растущие решения «сплющивают» базис. Чтобы этого не было, в некоторых
//! узлах базис переразлагается ([`Transform`]): каждый новый базисный вектор содержит ровно
//! одну растущую компоненту, частное решение — ни одной. Что такое «растущая компонента»,
//! определяет задача: вызывающий код передаёт матрицу растущих компонент базиса
//! (`r × cols`, коэффициенты разложения по растущим собственным векторам матрицы системы).
//! Цепочка преобразований ([`Chain`]) пересчитывает коэффициенты `c` назад, от правого конца к
//! левому.
//!
//! Базис хранится плоским построчным массивом `rows × cols`: строка — компонента вектора `y`,
//! столбец — решение; `cols = r` (однородная задача) или `r + 1` (последний столбец — частное
//! решение).

use crate::linalg::{self, CMatrix};
use num_complex::Complex64;

/// Когда переразлагать базис
#[derive(Debug, Clone, PartialEq)]
pub enum Policy {
    /// Никогда
    Never,
    /// Когда модуль действительной или мнимой части какого-либо элемента базиса превысил порог
    Threshold(f64),
    /// В заданных узлах, независимо от величины базиса (положения не зависят от параметров
    /// задачи — результат прогонки остаётся гладкой функцией этих параметров)
    Nodes(Vec<usize>),
}

impl Policy {
    /// Нужно ли переразложение в узле `node` при данном базисе
    pub fn due(&self, node: usize, basis: &[Complex64]) -> bool {
        match self {
            Policy::Never => false,
            Policy::Threshold(limit) => basis
                .iter()
                .any(|v| v.re.abs() > *limit || v.im.abs() > *limit),
            Policy::Nodes(nodes) => nodes.contains(&node),
        }
    }
}

/// Одно переразложение базиса: `X = Y·O`, `x₀ = y₀ − Y·o`
#[derive(Debug, Clone)]
pub struct Transform {
    /// `r × r`
    o: CMatrix,
    /// Длины `r`; для однородной задачи — нули
    shift: Vec<Complex64>,
}

impl Transform {
    /// Строит переразложение по растущим компонентам базиса
    /// # Arguments
    /// * `growing` - растущие компоненты базиса, построчно `r × cols`: элемент `(s, j)` —
    ///   коэффициент при `s`-й растущей компоненте в `j`-м столбце базиса
    /// * `r` - число базисных решений однородной системы
    /// * `cols` - число столбцов базиса: `r` или `r + 1` (с частным решением)
    /// * `target` - длины `r`: величина `s`-й растущей компоненты в новом `s`-м базисном
    ///   векторе (остальные растущие компоненты в нём равны нулю)
    /// # Result
    /// Returns `Ok(Transform)` с `O = H⁻¹·diag(target)`, `o = H⁻¹·h₀`, где `H` — первые `r`
    /// столбцов `growing`, `h₀` — последний (при `cols = r + 1`)
    /// # Errors
    /// Returns descriptive string if the sizes do not match or `H` is singular (the basis has
    /// lost linear independence in the growing components)
    pub fn new(
        growing: &[Complex64],
        r: usize,
        cols: usize,
        target: &[Complex64],
    ) -> Result<Transform, String> {
        check_cols(r, cols)?;
        if growing.len() != r * cols || target.len() != r {
            return Err(format!(
                "directed::Transform::new(): growing has {} elements (expected {}x{}), target has {} (expected {})",
                growing.len(),
                r,
                cols,
                target.len(),
                r
            ));
        }
        let h = CMatrix::from_fn(r, r, |s, j| growing[s * cols + j]);
        // Правые части: diag(target) и, если есть, растущие компоненты частного решения.
        let rhs = CMatrix::from_fn(r, cols, |s, j| {
            if j == r {
                growing[s * cols + r]
            } else if j == s {
                target[s]
            } else {
                Complex64::ZERO
            }
        });
        let solution = linalg::solve(&h, &rhs)
            .map_err(|e| format!("directed orthogonalization failed: {e}"))?;
        let shift = if cols > r {
            solution.column(r).iter().copied().collect()
        } else {
            vec![Complex64::ZERO; r]
        };
        Ok(Transform {
            o: solution.columns(0, r).into_owned(),
            shift,
        })
    }

    /// Число базисных решений
    pub fn rank(&self) -> usize {
        self.shift.len()
    }

    /// Преобразует базис на месте: `X = Y·O`, `x₀ = y₀ − Y·o`
    /// # Arguments
    /// * `basis` - базис, построчно `rows × cols`
    /// * `cols` - число столбцов: `r` или `r + 1`
    /// # Errors
    /// Returns descriptive string if the sizes do not match
    pub fn apply(&self, basis: &mut [Complex64], cols: usize) -> Result<(), String> {
        let r = self.rank();
        check_cols(r, cols)?;
        if cols == 0 || !basis.len().is_multiple_of(cols) {
            return Err(format!(
                "directed::Transform::apply(): basis of {} elements is not a matrix with {} columns",
                basis.len(),
                cols
            ));
        }
        let mut new_row = vec![Complex64::ZERO; cols];
        for row in basis.chunks_exact_mut(cols) {
            for (j, new) in new_row.iter_mut().enumerate().take(r) {
                *new = (0..r).map(|m| row[m] * self.o[(m, j)]).sum();
            }
            if cols > r {
                let shift: Complex64 = (0..r).map(|m| row[m] * self.shift[m]).sum();
                new_row[r] = row[r] - shift;
            }
            row.copy_from_slice(&new_row);
        }
        Ok(())
    }

    /// Коэффициенты разложения по базису **до** преобразования по коэффициентам **после**:
    /// `c_y = O·c_x − o`
    /// # Errors
    /// Returns descriptive string if `c` has a wrong length
    pub fn coefficients_before(&self, c: &[Complex64]) -> Result<Vec<Complex64>, String> {
        let r = self.rank();
        if c.len() != r {
            return Err(format!(
                "directed::Transform::coefficients_before(): {} coefficients, expected {}",
                c.len(),
                r
            ));
        }
        Ok((0..r)
            .map(|m| (0..r).map(|j| self.o[(m, j)] * c[j]).sum::<Complex64>() - self.shift[m])
            .collect())
    }
}

/// Цепочка переразложений вдоль прогонки
#[derive(Debug, Clone, Default)]
pub struct Chain {
    /// Узлы, в которых базис был преобразован (возрастают), и сами преобразования
    transforms: Vec<(usize, Transform)>,
}

impl Chain {
    pub fn new() -> Chain {
        Chain::default()
    }

    /// Запоминает преобразование, сделанное в узле `node`. Базис, сохранённый для этого узла,
    /// должен быть уже преобразованным.
    /// # Errors
    /// Returns descriptive string if `node` is not after the previous transform
    pub fn push(&mut self, node: usize, transform: Transform) -> Result<(), String> {
        if let Some((last, _)) = self.transforms.last()
            && node <= *last
        {
            return Err(format!(
                "directed::Chain::push(): node {node} is not after the previous transform at node {last}"
            ));
        }
        self.transforms.push((node, transform));
        Ok(())
    }

    /// Число переразложений
    pub fn len(&self) -> usize {
        self.transforms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.transforms.is_empty()
    }

    /// Коэффициенты для каждого участка прогонки
    /// # Arguments
    /// * `c_end` - коэффициенты разложения решения по базису на правом конце
    /// # Result
    /// Returns `Ok(segments)`: пары `(первый узел участка, коэффициенты)` по возрастанию узла,
    /// первый участок начинается в узле 0. Коэффициенты участка относятся к базисам всех узлов
    /// от его первого узла до первого узла следующего участка (не включая его).
    /// # Errors
    /// Returns descriptive string if `c_end` has a wrong length
    pub fn segments(&self, c_end: &[Complex64]) -> Result<Vec<(usize, Vec<Complex64>)>, String> {
        let mut segments = Vec::with_capacity(self.transforms.len() + 1);
        let mut c = c_end.to_vec();
        for (node, transform) in self.transforms.iter().rev() {
            let before = transform.coefficients_before(&c)?;
            segments.push((*node, std::mem::replace(&mut c, before)));
        }
        segments.push((0, c));
        segments.reverse();
        Ok(segments)
    }
}

/// Решение по базису и коэффициентам: `y = Y·c + y₀`
/// # Arguments
/// * `basis` - базис, построчно `rows × cols`
/// * `cols` - число столбцов: `c.len()` или `c.len() + 1` (с частным решением)
/// * `c` - коэффициенты
/// * `out` - решение, длины `rows`
/// # Errors
/// Returns descriptive string if the sizes do not match
pub fn combine(
    basis: &[Complex64],
    cols: usize,
    c: &[Complex64],
    out: &mut [Complex64],
) -> Result<(), String> {
    let r = c.len();
    check_cols(r, cols)?;
    if cols == 0 || basis.len() != out.len() * cols {
        return Err(format!(
            "directed::combine(): basis has {} elements, expected {}x{}",
            basis.len(),
            out.len(),
            cols
        ));
    }
    for (row, out) in basis.chunks_exact(cols).zip(out.iter_mut()) {
        let homogeneous: Complex64 = row.iter().zip(c).map(|(y, c)| y * c).sum();
        *out = if cols > r {
            homogeneous + row[r]
        } else {
            homogeneous
        };
    }
    Ok(())
}

/// Коэффициенты из условия на правом конце `D y = d`: `(D·Y)·c = d − D·y₀`
/// # Arguments
/// * `condition` - оператор условия, применённый к базису: `D·Y`, построчно `r × cols`
/// * `r` - число условий (и базисных решений)
/// * `cols` - число столбцов базиса: `r` или `r + 1`
/// * `d` - правая часть условия, длины `r`
/// # Result
/// Returns `Ok(c)`
/// # Errors
/// Returns descriptive string if the sizes do not match or the system is singular
pub fn solve_end_condition(
    condition: &[Complex64],
    r: usize,
    cols: usize,
    d: &[Complex64],
) -> Result<Vec<Complex64>, String> {
    check_cols(r, cols)?;
    if condition.len() != r * cols || d.len() != r {
        return Err(format!(
            "directed::solve_end_condition(): condition has {} elements (expected {}x{}), d has {} (expected {})",
            condition.len(),
            r,
            cols,
            d.len(),
            r
        ));
    }
    let matrix = CMatrix::from_fn(r, r, |s, j| condition[s * cols + j]);
    let rhs = CMatrix::from_fn(r, 1, |s, _| {
        if cols > r {
            d[s] - condition[s * cols + r]
        } else {
            d[s]
        }
    });
    let c = linalg::solve(&matrix, &rhs).map_err(|e| format!("end condition: {e}"))?;
    Ok(c.iter().copied().collect())
}

fn check_cols(r: usize, cols: usize) -> Result<(), String> {
    if cols == r || cols == r + 1 {
        Ok(())
    } else {
        Err(format!(
            "directed: basis of {r} solutions must have {r} or {} columns, got {cols}",
            r + 1
        ))
    }
}

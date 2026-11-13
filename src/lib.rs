pub mod algebra;
pub mod bvp;
pub mod gsl;
pub mod linalg;
pub mod mesh;
pub mod ode;

/// Плотные матрицы и разложения: типы `nalgebra` встречаются в API (`linalg`).
pub use nalgebra;
/// Комплексные числа, общие для всех пользователей библиотеки.
pub use num_complex;

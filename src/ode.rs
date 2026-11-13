//! Задача Коши для систем обыкновенных дифференциальных уравнений.

use std::ops::{Add, Mul};

/// Classical fourth-order Runge–Kutta stepper with preallocated work buffers.
///
/// The state is a flat slice of `T` (`f64`, `Complex64`, …); a matrix of several solutions is
/// integrated as one flat state.
pub struct Rk4<T> {
    k1: Vec<T>,
    k2: Vec<T>,
    k3: Vec<T>,
    k4: Vec<T>,
    tmp: Vec<T>,
}

impl<T> Rk4<T>
where
    T: Copy + Default + Add<Output = T> + Mul<f64, Output = T>,
{
    /// Creates a stepper for a state of `n` elements
    pub fn new(n: usize) -> Rk4<T> {
        let buffer = || vec![T::default(); n];
        Rk4 {
            k1: buffer(),
            k2: buffer(),
            k3: buffer(),
            k4: buffer(),
            tmp: buffer(),
        }
    }

    /// One step from `x0` to `x1` (the step may be negative)
    /// # Arguments
    /// * `rhs` - right-hand side `rhs(x, y, dy)`: writes `dy/dx` into `dy` (every element)
    /// * `x0`, `x1` - ends of the step
    /// * `y` - state at `x0` on entry, at `x1` on exit
    /// # Errors
    /// Returns descriptive string if the length of `y` differs from the one given to `new`,
    /// or the error of `rhs` (the state is then left unchanged)
    pub fn step(
        &mut self,
        rhs: &mut impl FnMut(f64, &[T], &mut [T]) -> Result<(), String>,
        x0: f64,
        x1: f64,
        y: &mut [T],
    ) -> Result<(), String> {
        if y.len() != self.tmp.len() {
            return Err(format!(
                "Rk4::step(): state has {} elements, the stepper was created for {}",
                y.len(),
                self.tmp.len()
            ));
        }
        let h = x1 - x0;
        let xm = x0 + 0.5 * h;
        rhs(x0, y, &mut self.k1)?;
        advance(&mut self.tmp, y, &self.k1, 0.5 * h);
        rhs(xm, &self.tmp, &mut self.k2)?;
        advance(&mut self.tmp, y, &self.k2, 0.5 * h);
        rhs(xm, &self.tmp, &mut self.k3)?;
        advance(&mut self.tmp, y, &self.k3, h);
        rhs(x1, &self.tmp, &mut self.k4)?;
        let (k1, k2, k3, k4) = (&self.k1, &self.k2, &self.k3, &self.k4);
        for (i, y) in y.iter_mut().enumerate() {
            *y = *y + (k1[i] + (k2[i] + k3[i]) * 2.0 + k4[i]) * (h / 6.0);
        }
        Ok(())
    }
}

/// `out = y + k * h`
fn advance<T>(out: &mut [T], y: &[T], k: &[T], h: f64)
where
    T: Copy + Add<Output = T> + Mul<f64, Output = T>,
{
    for ((out, &y), &k) in out.iter_mut().zip(y).zip(k) {
        *out = y + k * h;
    }
}

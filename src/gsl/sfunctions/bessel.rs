mod gsl_signatures;

pub fn jn(n: i32, x: f64) -> Result<f64, String> {
    Ok(unsafe { gsl_signatures::gsl_sf_bessel_Jn(n, x) })
}

pub fn djn(n: i32, x: f64) -> Result<f64, String> {
    Ok((jn(n - 1, x)? - jn(n + 1, x)?) * 0.5)
}

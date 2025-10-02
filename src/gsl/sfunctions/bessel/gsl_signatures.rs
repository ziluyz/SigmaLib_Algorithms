use std::os::raw::{c_double, c_int};

unsafe extern "C" {
    pub fn gsl_sf_bessel_Jn(n: c_int, x: c_double) -> c_double;
}

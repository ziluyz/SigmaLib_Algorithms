use std::os::raw::{c_double, c_int};

/// `gsl_sf_result`
#[repr(C)]
pub struct GslSfResult {
    pub val: c_double,
    pub err: c_double,
}

unsafe extern "C" {
    pub fn gsl_sf_bessel_Jn_e(n: c_int, x: c_double, result: *mut GslSfResult) -> c_int;
}

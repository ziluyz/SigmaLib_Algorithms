use super::GslFunction;

unsafe extern "C" {
    pub fn gsl_integration_qng(
        f: *const GslFunction,
        a: f64,
        b: f64,
        epsabs: f64,
        epsrel: f64,
        result: *mut f64,
        abserr: *mut f64,
        neval: *mut usize,
    ) -> i32;
}

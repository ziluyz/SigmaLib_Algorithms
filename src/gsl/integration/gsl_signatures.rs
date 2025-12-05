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

    pub fn gsl_integration_workspace_alloc(n: usize) -> *mut std::ffi::c_void;
    pub fn gsl_integration_workspace_free(workspace: *mut std::ffi::c_void);

    pub fn gsl_integration_qawo_table_alloc(
        omega: f64,
        L: f64,
        sine: i32,
        n: usize,
    ) -> *mut std::ffi::c_void;
    pub fn gsl_integration_qawo_table_set(
        t: *mut std::ffi::c_void,
        omega: f64,
        L: f64,
        sine: i32,
    ) -> i32;
    pub fn gsl_integration_qawo_table_set_length(t: *mut std::ffi::c_void, L: f64) -> i32;
    pub fn gsl_integration_qawo_table_free(t: *mut std::ffi::c_void);
    pub fn gsl_integration_qawo(
        f: *const GslFunction,
        a: f64,
        epsabs: f64,
        epsrel: f64,
        limit: usize,
        workspace: *mut std::ffi::c_void,
        wf: *mut std::ffi::c_void,
        result: *mut f64,
        abserr: *mut f64,
    ) -> i32;

    pub fn gsl_integration_romberg_alloc(n: usize) -> *mut std::ffi::c_void;
    pub fn gsl_integration_romberg_free(workspace: *mut std::ffi::c_void);
    pub fn gsl_integration_romberg(
        f: *const GslFunction,
        a: f64,
        b: f64,
        epsabs: f64,
        epsrel: f64,
        result: *mut f64,
        neval: *mut usize,
        workspace: *mut std::ffi::c_void,
    ) -> i32;
}

mod gsl_signatures;
mod qng;

#[repr(C)]
pub struct GslFunction {
    pub function: Option<
        unsafe extern "C" fn(
            x: std::os::raw::c_double,
            params: *mut std::ffi::c_void,
        ) -> std::os::raw::c_double,
    >,
    pub params: *mut std::ffi::c_void,
}

unsafe extern "C" fn trampoline(x: f64, params: *mut std::ffi::c_void) -> f64 {
    let closure = unsafe { &mut *(params as *mut &mut dyn FnMut(f64) -> f64) };
    closure(x)
}

pub use qng::{QngResult, qng};

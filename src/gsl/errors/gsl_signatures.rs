unsafe extern "C" {
    pub fn gsl_set_error_handler_off() -> *const std::ffi::c_void;
    pub fn gsl_strerror(status: std::os::raw::c_int) -> *const std::os::raw::c_char;
}

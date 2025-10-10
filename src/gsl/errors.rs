mod gsl_signatures;

pub fn turn_off_error_handler() {
    unsafe { gsl_signatures::gsl_set_error_handler_off() };
}

pub fn get_error_description(status: i32) -> String {
    unsafe {
        let err_ptr = gsl_signatures::gsl_strerror(status);
        let c_str = std::ffi::CStr::from_ptr(err_ptr);
        c_str.to_string_lossy().into_owned()
    }
}

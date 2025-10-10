use super::{GslFunction, gsl_signatures, trampoline};

pub struct QngResult {
    pub result: f64,
    pub abserr: f64,
    pub neval: usize,
}

pub fn qng<F>(f: &mut F, a: f64, b: f64, epsabs: f64, epsrel: f64) -> Result<QngResult, String>
where
    F: FnMut(f64) -> f64,
{
    let mut result = 0.0;
    let mut abserr = 0.0;
    let mut neval: usize = 0;

    // Приводим замыкание к типу, который понимает trampoline
    let mut closure_trait_obj: &mut dyn FnMut(f64) -> f64 = f;
    let closure_ptr = &mut closure_trait_obj as *mut _ as *mut std::ffi::c_void;
    let gsl_func = GslFunction {
        function: Some(trampoline),
        params: closure_ptr,
    };

    let status = unsafe {
        gsl_signatures::gsl_integration_qng(
            &gsl_func,
            a,
            b,
            epsabs,
            epsrel,
            &mut result,
            &mut abserr,
            &mut neval,
        )
    };

    if status == 0 {
        Ok(QngResult {
            result,
            abserr,
            neval,
        })
    } else {
        Err(format!(
            "Error in gsl_integration_qng: {}",
            super::super::errors::get_error_description(status)
        ))
    }
}

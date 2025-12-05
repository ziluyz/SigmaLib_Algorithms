use super::{GslFunction, gsl_signatures, trampoline};

pub struct RombergResult {
    pub result: f64,
    pub neval: usize,
}

pub struct Romberg {
    workspace: *mut std::ffi::c_void,
    a: f64,
    b: f64,
    iterations: Option<usize>,
}

impl Romberg {
    pub fn new(a: f64, b: f64) -> Romberg {
        Romberg {
            workspace: std::ptr::null_mut(),
            a,
            b,
            iterations: None,
        }
    }

    pub fn initialize(&mut self, iterations: usize) {
        unsafe {
            self.workspace = gsl_signatures::gsl_integration_romberg_alloc(iterations);
            self.iterations = Some(iterations);
        }
    }

    pub fn integrate<F>(&self, f: &mut F, epsabs: f64, epsrel: f64) -> Result<RombergResult, String>
    where
        F: FnMut(f64) -> f64,
    {
        if self.iterations.is_none() {
            return Err("Romberg is not initialized".to_string());
        }

        let mut result = 0.0;
        let mut neval = 0;

        let mut closure_trait_obj: &mut dyn FnMut(f64) -> f64 = f;
        let closure_ptr = &mut closure_trait_obj as *mut _ as *mut std::ffi::c_void;
        let gsl_func = GslFunction {
            function: Some(trampoline),
            params: closure_ptr,
        };

        let status = unsafe {
            gsl_signatures::gsl_integration_romberg(
                &gsl_func,
                self.a,
                self.b,
                epsabs,
                epsrel,
                &mut result,
                &mut neval,
                self.workspace,
            )
        };

        if status == 0 {
            Ok(RombergResult { result, neval })
        } else {
            Err(format!(
                "Error in gsl_integration_romberg: {}",
                super::super::errors::get_error_description(status)
            ))
        }
    }
}

impl Drop for Romberg {
    fn drop(&mut self) {
        if self.iterations.is_none() {
            return;
        }
        unsafe {
            gsl_signatures::gsl_integration_romberg_free(self.workspace);
        }
    }
}

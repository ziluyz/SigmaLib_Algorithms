use super::{GslFunction, gsl_signatures, trampoline};

pub struct QawoResult {
    pub result: f64,
    pub abserr: f64,
}

pub enum QawoOscillator {
    Sine,
    Cosine,
}

impl QawoOscillator {
    fn interpret(&self) -> i32 {
        match self {
            Self::Sine => 1,
            Self::Cosine => 0,
        }
    }
}

pub struct Qawo {
    workspace: *mut std::ffi::c_void,
    table: *mut std::ffi::c_void,
    a: f64,
    b: f64,
    intervals: Option<usize>,
}

impl Qawo {
    pub fn new(a: f64, b: f64) -> Qawo {
        Qawo {
            workspace: std::ptr::null_mut(),
            table: std::ptr::null_mut(),
            a,
            b,
            intervals: None,
        }
    }

    pub fn initialize(&mut self, omega: f64, intervals: usize, sine: QawoOscillator) {
        unsafe {
            self.workspace = gsl_signatures::gsl_integration_workspace_alloc(intervals);
            let log2 = ((intervals as f64).log(2.0) as usize).max(2);
            self.table = gsl_signatures::gsl_integration_qawo_table_alloc(
                omega,
                self.b - self.a,
                sine.interpret(),
                log2,
            );
            self.intervals = Some(intervals);
        }
    }

    pub fn change_parameters(&mut self, omega: f64, a: f64, b: f64, sine: QawoOscillator) {
        self.a = a;
        self.b = b;
        unsafe {
            gsl_signatures::gsl_integration_qawo_table_set(
                self.table,
                omega,
                b - a,
                sine.interpret(),
            );
        }
    }

    pub fn change_limits(&mut self, a: f64, b: f64) {
        self.a = a;
        self.b = b;
        unsafe {
            gsl_signatures::gsl_integration_qawo_table_set_length(self.table, b - a);
        }
    }

    pub fn integrate<F>(
        &self,
        f: &mut F,
        epsabs: f64,
        epsrel: f64,
        ignore_roundoff_error: bool,
    ) -> Result<QawoResult, String>
    where
        F: FnMut(f64) -> f64,
    {
        if self.intervals.is_none() {
            return Err("QAWO is not initialized".to_string());
        }

        let mut result = 0.0;
        let mut abserr = 0.0;

        let mut closure_trait_obj: &mut dyn FnMut(f64) -> f64 = f;
        let closure_ptr = &mut closure_trait_obj as *mut _ as *mut std::ffi::c_void;
        let gsl_func = GslFunction {
            function: Some(trampoline),
            params: closure_ptr,
        };

        let status = unsafe {
            gsl_signatures::gsl_integration_qawo(
                &gsl_func,
                self.a,
                epsabs,
                epsrel,
                self.intervals.unwrap(),
                self.workspace,
                self.table,
                &mut result,
                &mut abserr,
            )
        };

        if status == 0 || (status == 18 && ignore_roundoff_error) {
            Ok(QawoResult { result, abserr })
        } else {
            Err(format!(
                "Error in gsl_integration_qawo: {}",
                super::super::errors::get_error_description(status)
            ))
        }
    }
}

impl Drop for Qawo {
    fn drop(&mut self) {
        if self.intervals.is_none() {
            return;
        }
        unsafe {
            gsl_signatures::gsl_integration_workspace_free(self.workspace);
            gsl_signatures::gsl_integration_qawo_table_free(self.table);
        }
    }
}

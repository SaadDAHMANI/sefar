use crate::core::OptError;

pub trait Parameters {
    fn get_population_size(&self) -> usize {
        10usize
    }

    fn get_problem_dimension(&self) -> usize {
        10usize
    }

    /// The defaut stopping citerion is Max-Iteration = 1.
    fn get_max_iterations(&self) -> usize {
        1usize
    }

    fn get_lower_bounds(&self) -> Vec<f64> {
        let d = self.get_problem_dimension();
        let mut lb = Vec::new();
        for _i in 0..d {
            lb.push(-100.0f64);
        }
        lb
    }

    fn get_upper_bounds(&self) -> Vec<f64> {
        let d = self.get_problem_dimension();
        let mut ub = Vec::new();
        for _i in 0..d {
            ub.push(100.0f64);
        }
        ub
    }

    #[cfg(feature = "parallel")]
    fn get_num_threads(&self) -> ThreadPoolSize {
        ThreadPoolSize::Auto
    }

    fn check(&self) -> Result<(), OptError> {
        if self.get_population_size() == 0 {
            return Err(OptError::PopulationSizeIsNull);
        }

        if self.get_problem_dimension() == 0 {
            return Err(OptError::ProblemDimensionIsNull);
        }

        if self.get_max_iterations() == 0 {
            return Err(OptError::MaxIterationsIsNull);
        }

        if self.get_lower_bounds().is_empty() {
            return Err(OptError::EmptyLB);
        }

        if self.get_upper_bounds().is_empty() {
            return Err(OptError::EmptyUB);
        }

        if self.get_lower_bounds().len() != self.get_upper_bounds().len() {
            return Err(OptError::LBLengthNotEqualsUBlength {
                lb_len: self.get_lower_bounds().len(),
                ub_len: self.get_upper_bounds().len(),
            });
        }

        if self.get_lower_bounds().len() != self.get_problem_dimension() {
            return Err(OptError::LBLengthNoEqualsProblemDim {
                lb_len: self.get_lower_bounds().len(),
                problem_dim: self.get_problem_dimension(),
            });
        }

        if self.get_upper_bounds().len() != self.get_problem_dimension() {
            return Err(OptError::UBLengthNoEqualsProblemDim {
                ub_len: self.get_lower_bounds().len(),
                problem_dim: self.get_problem_dimension(),
            });
        }

        Ok(())
    }
}

/// Use this enum to define number of worker threads for parallel run.
/// Auto : use the maximum number of threads (all available cores)
/// Fixed(n): use 'n' threads.  
/// # Example
///
/// ```
/// let num_threads = ThreadPoolSize::Auto;
/// // or
/// let num_threads = ThreadPoolSize::Fixed(2);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ThreadPoolSize {
    Auto,
    Fixed(usize),
}

impl ThreadPoolSize {
    pub fn as_usize(&self) -> usize {
        match self {
            Self::Auto => rayon::current_num_threads(),
            Self::Fixed(n) => {
                let nb_thr = *n;
                nb_thr.min(rayon::current_num_threads())
            }
        }
    }
}

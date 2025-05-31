use crate::core::OptError;

pub trait Parameters {
    fn get_population_size(&self) -> usize {
        10usize
    }

    fn get_problem_dimension(&self) -> usize {
        10usize
    }

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
            return Err(OptError::LBLengthNotEqualsUBlength);
        }

        if self.get_lower_bounds().len() != self.get_problem_dimension() {
            return Err(OptError::LBLengthNoEqualsProblemDim);
        }

        if self.get_upper_bounds().len() != self.get_problem_dimension() {
            return Err(OptError::LBLengthNoEqualsProblemDim);
        }

        Ok(())
    }
}

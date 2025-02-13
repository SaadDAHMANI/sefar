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

    fn check(&self) -> Result<(), String> {
        let mut errors: usize = 0;
        let mut msg: String = String::new();

        if self.get_population_size() == 0 {
            msg = String::from("population_size must be greater than 0!; \n");
            errors += 1;
        }

        if self.get_problem_dimension() == 0 {
            msg = format!(
                "{} Search space dimension (i.e., problem dimension or decision variables) must be greater than 0!; \n",
                msg
            );
            errors += 1;
        }

        if self.get_max_iterations() == 0 {
            msg = format!(
                "{} Iterations count (max_iterations) must be greater than 0!; \n",
                msg
            );
            errors += 1;
        }

        if self.get_lower_bounds().is_empty() {
            msg = format!("{} Lower_bounds length must be greater than 0!; \n", msg);
            errors += 1;
        }

        if self.get_upper_bounds().is_empty() {
            msg = format!("{} Upper_bounds length must be greater than 0!; \n", msg);
            errors += 1;
        }

        if self.get_lower_bounds().len() != self.get_upper_bounds().len() {
            msg = format!(
                "{} Lower_bounds & Upper_bounds lengths must be equal!; \n",
                msg
            );
            errors += 1;
        }

        if self.get_lower_bounds().len() != self.get_problem_dimension()
            || self.get_upper_bounds().len() != self.get_problem_dimension()
        {
            msg = format!(
                "{} Lower_bounds & Upper_bounds lengths must equal search space dimension!; \n",
                msg
            );
            errors += 1;
        }

        if errors > 0 {
            msg = format!("There are [{}] errors : \n {}", errors, msg.trim());
            Err(msg)
        } else {
            Ok(())
        }
    }
}

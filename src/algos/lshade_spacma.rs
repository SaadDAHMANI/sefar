use crate::common::*;
use crate::core::eoa::{InitializationMode, EOA};
use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
//use rand::rngs::ThreadRng;
use rand_distr::{Distribution, Uniform};
//#[cfg(feature = "parallel")]
//use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use std::time::Instant;

/// LshadeSpacma(LSHADE_SPACMA)
/// Reference:
/// Ali W. Mohamed, Anas A. Hadi, Anas M. Fattouh, and Kamal M. Jambi:
/// L-SHADE with Semi Parameter Adaptation Approach for Solving CEC 2017 Benchmark Problems,
/// Proc. IEEE Congress on Evolutionary Computation (CEC-2017), Spain, June, 2017
/// https://ieeexplore.ieee.org/document/7969307

#[derive(Debug)]
pub struct LshadeSpacma<'a, T: Problem> {
    /// The problem to optimize. It must define the Problem trait.
    pub problem: &'a mut T,

    /// Define the parameters of GO algorithm.
    pub params: &'a LshadeSpacmaParams<'a>,
}

impl<'a, T: Problem> EOA for LshadeSpacma<'a, T> {
    fn run(&mut self) -> OptimizationResult {
        let result: OptimizationResult = OptimizationResult::get_none(String::from("n/a"));
        let l_rate: f64 = 0.80;
        let record_fes_factor: Vec<f64> = vec![
            0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        ];

        let progress: usize = record_fes_factor.len();
        //-----------------------------------------------------------------
        let problem_size: usize = self.params.get_problem_dimension();
        let pop_size: usize = self.params.get_population_size();
        let max_iter: usize = self.params.get_max_iterations();
        let max_nfes: usize = (max_iter + 1) * pop_size;
        let mut nfes: usize = 0;
        let mut run_funcvals: Vec<f64> = vec![-1.0; max_iter + 1];

        //
        // Parameter settings for L-SHADE----------------------------------
        let p_best_rate: f64 = 0.11;
        let arc_rate: f64 = 1.4;
        let memory_size: usize = 5;
        // pop_size = 18 * problem_size;
        let max_pop_size: usize = usize::max(pop_size, 18 * problem_size);
        let min_pop_size: usize = 4;
        //------------------------------------------------------------------
        result
    }
}

#[derive(Debug, Clone)]
pub struct LshadeSpacmaParams<'a> {
    /// The number of search agents.
    pub population_size: usize,

    /// The dimension of the optimization problem (i.e., the length of the solution).
    pub dimensions: usize,

    /// The maximum number of iterations serves as the stopping criterion for the optimization process.
    pub max_iterations: usize,

    /// The lower bounds of the search space.
    pub lower_bounds: &'a [f64],

    /// The upper bounds of the search space.
    pub upper_bounds: &'a [f64],
}

impl<'a> Parameters for LshadeSpacmaParams<'a> {
    fn get_problem_dimension(&self) -> usize {
        self.dimensions
    }

    fn get_max_iterations(&self) -> usize {
        usize::max(self.max_iterations, 1)
    }

    fn get_population_size(&self) -> usize {
        usize::max(self.population_size, 12)
    }

    fn get_lower_bounds(&self) -> Vec<f64> {
        self.lower_bounds.to_vec()
    }

    fn get_upper_bounds(&self) -> Vec<f64> {
        self.upper_bounds.to_vec()
    }
}

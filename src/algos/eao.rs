// EAO : Enzyme action optimizer

use crate::common::*;
use crate::core::eoa::{InitializationMode, EOA};
use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
use rand::rngs::ThreadRng;
use rand_distr::{Distribution, Uniform};
// #[cfg(feature = "parallel")]
// use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use std::fmt::Display;
use std::time::Instant;

///
/// EAO : Enzyme action optimizer
///
/// Reference:
/// Enzyme action optimizer: a novel bio-inspired optimization algorithm
/// Rodan, A., Al-Tamimi, A. K., Al-Alnemer, L., Mirjalili, S., & Tiňo, P. (2025).
/// Enzyme action optimizer: a novel bio-inspired optimization algorithm.
/// The Journal of Supercomputing, 81(5), 686
/// Paper link: <https://link.springer.com/article/10.1007/s11227-025-07052-w>
/// Original source code : <https://github.com/AliRodan/Enzyme-Action-Optimizer>
///  
/// Written in Rust by Saad Dahmani <sd.dahmani2000@gmail.com>
///

///
/// Define parameters for the Growth (GO) algorithm.
///
#[derive(Debug, Clone)]
pub struct EAOparams<'a> {
    /// The number of search agents.
    pub population_size: usize,

    /// The dimension of the optimization problem (i.e., the length of the solution).
    pub problem_dimension: usize,

    /// The maximum number of iterations serves as the stopping criterion for the optimization process.
    pub max_iterations: usize,

    /// The lower bounds of the search space.
    pub lower_bounds: &'a [f64],

    /// The upper bounds of the search space.
    pub upper_bounds: &'a [f64],
}

impl<'a> EAOparams<'a> {
    ///
    /// Create a new instance of GO parameters:
    /// pop_size : The number of search agents.
    /// dim : The dimension of the optimization problem.
    /// max_iter : The maximum number of iterations serves as the stopping criterion.
    /// lb : The lower bounds of the search space.
    /// ub : The upper bounds of the search space.
    ///
    #[allow(dead_code)]
    pub fn new(pop_size: usize, dim: usize, max_iter: usize, lb: &'a [f64], ub: &'a [f64]) -> Self {
        EAOparams {
            population_size: pop_size,
            problem_dimension: dim,
            max_iterations: max_iter,
            lower_bounds: lb,
            upper_bounds: ub,
        }
    }
}

impl<'a> Parameters for EAOparams<'a> {
    fn get_problem_dimension(&self) -> usize {
        self.problem_dimension
    }

    fn get_max_iterations(&self) -> usize {
        usize::max(self.max_iterations, 1)
    }

    fn get_population_size(&self) -> usize {
        usize::max(self.population_size, 6)
    }

    fn get_lower_bounds(&self) -> Vec<f64> {
        self.lower_bounds.to_vec()
    }

    fn get_upper_bounds(&self) -> Vec<f64> {
        self.upper_bounds.to_vec()
    }
}

impl<'a> Default for EAOparams<'a> {
    ///
    /// Return the default values of parameters, as follows:
    ///
    /// ~~~
    ///
    ///  use sefar::algos::go::*;
    ///
    ///  GOparams {
    ///     population_size : 20,
    ///     problem_dimension : 3,
    ///     max_iterations : 100,
    ///     lower_bounds : &[100.0f64, 100.0, 100.0],
    ///     upper_bounds : &[-100.0f64, -100.0, -100.0],
    /// };
    /// ~~~
    ///
    fn default() -> Self {
        EAOparams {
            population_size: 20,
            problem_dimension: 3,
            max_iterations: 100,
            lower_bounds: &[-100.0f64, -100.0, -100.0],
            upper_bounds: &[100.0f64, 100.0, 100.0],
        }
    }
}

impl<'a> Display for EAOparams<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Popo.Size: {}, Problem dim.: {}, Max.Iter: {}, LB: {:?}, UB: {:?}",
            self.population_size,
            self.problem_dimension,
            self.max_iterations,
            self.get_lower_bounds(),
            self.get_upper_bounds()
        )
    }
}

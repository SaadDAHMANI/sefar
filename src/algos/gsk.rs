use crate::common::*;
use crate::core::eoa::{InitializationMode, EOA};
use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
use rand::rngs::ThreadRng;
use rand_distr::{Distribution, Uniform};
//#[cfg(feature = "parallel")]
//use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use std::time::Instant;

/// GSK : Gaining-Sharing Knowedge algorithm.
/// Reference:
/// Mohamed, A. W., Hadi, A. A., & Mohamed, A. K. (2020).
/// Gaining-sharing knowledge based algorithm for solving optimization problems: a novel nature-inspired algorithm.
/// International Journal of Machine Learning and Cybernetics, 11(7), 1501-1529.
/// (https​://doi.org/10.1007/s1304​2-019-01053​-x)
/// Matlab original code : <https://sites.google.com/view/optimization-project/files?authuser=0>
///
/// Written in Rust by Saad Dahmani <sd.dahmani2000@gmail.com>
///

#[derive(Debug)]
pub struct GSK<'a, T: Problem> {
    /// The problem to optimize. It must define the Problem trait.
    pub problem: &'a mut T,

    /// Define the parameters of GO algorithm.
    pub params: &'a GSKparams<'a>,
}

impl<'a, T: Problem> GSK<'a, T> {
    ///
    /// Return a new instance of the Gaining-Sharing Knowledge Optimizer (GSK).
    /// settings: The optimization parameters,
    /// problem: The problem to optimize.
    ///
    pub fn new(settings: &'a GSKparams, problem: &'a mut T) -> Self {
        GSK {
            problem,
            params: settings,
        }
    }
}

impl<'a, T: Problem> EOA for GSK<'a, T> {
    fn run(&mut self) -> OptimizationResult {
        let result: OptimizationResult = OptimizationResult::get_none(String::from("n/a"));

        //-------------------------------------------------
        let pop_size: usize = self.params.population_size;
        let max_iter: usize = self.params.max_iterations;
        let dim: usize = self.params.dimensions;
        let max_nfes: usize = pop_size * (max_iter + 1);
        //--------------------------------------------------
        let mut nfes: usize = 0; // function evaluation counter.
        let mut bsf_fit_var: f64 = f64::MAX; // the best fitness value.
        let mut fitness: Vec<f64> = vec![0.0f64; pop_size];
        let mut run_funcvals: Vec<f64> = vec![0.0f64; max_iter];
        //--------------------------------------------------

        let g_max_f64: f64 = max_nfes as f64 / pop_size as f64;
        let g_max: usize = g_max_f64.trunc() as usize;

        // Initialize the main population:
        // Initialize the old population
        let mut popold = self.initialize(self.params, InitializationMode::RealUniform);

        // Initialize the current population
        let mut pop = popold.clone();

        // Objective function evaluation:
        for i in 0..pop_size {
            fitness[i] = self.problem.objectivefunction(&pop[i].genes);
            nfes += 1;
            //println!("fitness[{}] = {}", i, fitness[i]);
        }
        // Save the best fitness value for convergence trend:
        for i in 0..pop_size {
            if fitness[i] < bsf_fit_var {
                bsf_fit_var = fitness[i];
            }
        }
        run_funcvals[0] = bsf_fit_var; //save history of convergence.

        //--------------------------------------------------
        let kf = 0.5; //Knowledge Factor.
        let kr = 0.9; //Knowledge Ratio.
        let k = vec![10.0; pop_size];
        let g: usize = 0;

        //THE MAIN LOOP

        //--------------------------------------------------

        result
    }
}

///
/// Define parameters for the Gaining-Sharing Knowledge (GSK) Algorithm.
///
#[derive(Debug, Clone)]
pub struct GSKparams<'a> {
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

impl<'a> Parameters for GSKparams<'a> {
    fn get_dimensions(&self) -> usize {
        self.dimensions
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

impl<'a> Default for GSKparams<'a> {
    ///
    /// Return the default values of parameters, as follows:
    ///
    /// ~~~
    ///
    ///  use sefar::algos::gsk::*;
    ///
    ///  GOparams {
    ///     population_size : 12,
    ///     dimensions : 3,
    ///     max_iterations : 1,
    ///     lower_bounds : &[100.0f64, 100.0, 100.0],
    ///     upper_bounds : &[-100.0f64, -100.0, -100.0],
    /// };
    /// ~~~
    ///
    fn default() -> Self {
        GSKparams {
            population_size: 12,
            dimensions: 3,
            max_iterations: 1,
            lower_bounds: &[-100.0f64, -100.0, -100.0],
            upper_bounds: &[100.0f64, 100.0, 100.0],
        }
    }
}

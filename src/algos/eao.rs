// EAO : Enzyme action optimizer

use crate::common::*;
use crate::core::eoa::{InitializationMode, EOA};
use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
use rand::rngs::ThreadRng;
use rand_distr::num_traits::real::Real;
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
#[derive(Debug)]
pub struct EAO<'a, T: Problem> {
    /// The problem to optimize. It must define the Problem trait.
    pub problem: &'a mut T,

    /// Define the parameters of EAO algorithm.
    pub params: &'a EAOparams<'a>,
}

impl<'a, T: Problem> EAO<'a, T> {
    ///
    /// Return a new instance of Enzyme Action Optimizer.
    /// settings: Optimization parameters,
    /// problem: Problem to optimize.
    ///
    pub fn new(settings: &'a EAOparams, problem: &'a mut T) -> Self {
        EAO {
            problem,
            params: settings,
        }
    }
}

impl<'a, T: Problem> EOA for EAO<'a, T> {
    ///
    /// Call this function to execute GO algorithm.
    ///
    fn run(&mut self) -> OptimizationResult {
        // let chronos = Instant::now();

        let enzyme_count: usize = self.params.population_size;
        let active_site_dimension: usize = self.params.problem_dimension;

        let max_iter: usize = self.params.max_iterations;
        let mut break_process: bool = false;

        let ub = self.params.upper_bounds;
        let lb = self.params.lower_bounds;

        let mut reaction_rate: Vec<f64> = vec![0.0; enzyme_count];

        let mut convergence_curve = vec![0.0f64; max_iter];

        let mut af: f64 = 0.0;
        // =============================================================
        let mut substrate_pool = self.initialize(self.params, InitializationMode::RealUniform);

        // ______________ fitness evaluation ___________________________
        //Evaluation of search agents
        // Sequential mode
        #[cfg(not(feature = "parallel"))]
        for i in 0..enzyme_count {
            reaction_rate[i] = self.problem.objectivefunction(&mut substrate_pool[i].genes);
        }

        //___________Parallel mode________________
        #[cfg(feature = "parallel")]
        {
            substrate_pool
                .par_iter_mut()
                .for_each(|g| g.fitness = Some(self.problem.objectivefunction(&g.genes)));
            for i in 0..n {
                match substrate_pool[i].fitness {
                    None => substrate_pool[i] = f64::MAX,
                    Some(fit) => substrate_pool[i] = fit,
                };
            }
        }
        //________________________________________

        // [OptimalCatalysis, idx] = min(ReactionRate);

        let mut best_substrate: f64 = reaction_rate
            .iter()
            .copied()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(f64::MAX);
        // _____________________________________________________________
        // =============================================================

        // MAIN LOOP
        for t in 1..max_iter + 1 {
            af = f64::sqrt(t as f64 / max_iter as f64);

            for i in 0..enzyme_count {
                //  1) Update FirstSubstratePosition
            }
        }

        // Return result

        OptimizationResult::get_empty(None)
    }
}

///
/// Define parameters for the Growth (EAO) algorithm.
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

    /// Enzyme Concentration EC. The default value = 0.1
    pub ec: f64,
}

impl<'a> EAOparams<'a> {
    ///
    /// Create a new instance of EAO parameters:
    /// pop_size : The number of search agents.
    /// dim : The dimension of the optimization problem.
    /// max_iter : The maximum number of iterations serves as the stopping criterion.
    /// lb : The lower bounds of the search space.
    /// ub : The upper bounds of the search space.
    /// ec : Enzyme concentration, default value = 0.1.
    ///
    #[allow(dead_code)]
    pub fn new(
        pop_size: usize,
        dim: usize,
        max_iter: usize,
        lb: &'a [f64],
        ub: &'a [f64],
        ec: f64,
    ) -> Self {
        EAOparams {
            population_size: pop_size,
            problem_dimension: dim,
            max_iterations: max_iter,
            lower_bounds: lb,
            upper_bounds: ub,
            ec,
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
    ///  EAOparams {
    ///     population_size : 20,
    ///     problem_dimension : 3,
    ///     max_iterations : 100,
    ///     lower_bounds : &[100.0f64, 100.0, 100.0],
    ///     upper_bounds : &[-100.0f64, -100.0, -100.0],
    /// ec : 0.1,
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
            ec: 0.1,
        }
    }
}

impl<'a> Display for EAOparams<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Popo.Size: {}, Problem dim.: {}, Max.Iter: {}, LB: {:?}, UB: {:?}, EC: {:?}",
            self.population_size,
            self.problem_dimension,
            self.max_iterations,
            self.get_lower_bounds(),
            self.get_upper_bounds(),
            self.ec
        )
    }
}

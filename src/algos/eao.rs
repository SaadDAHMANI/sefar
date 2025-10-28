// EAO : Enzyme action optimizer

use crate::common::*;
use crate::core::eoa::{InitializationMode, EOA};
use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
// use rand::rngs::ThreadRng;
// use rand_distr::num_traits::real::Real;
use rand_distr::{Distribution, Uniform};
#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
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

    fn pick_two_random_distinct_substrates(&self) -> (usize, usize) {
        let between = Uniform::from(0..self.params.population_size);
        let mut rng = rand::thread_rng();

        let indx1: usize = between.sample(&mut rng);
        let mut indx2: usize = indx1;
        let mut counter: usize = 0;

        while indx1 == indx2 {
            indx2 = between.sample(&mut rng);
            counter += 1;
            if counter > 100 {
                break;
            }
            // println!("counter = {}", counter);
        }
        (indx1, indx2)
    }

    fn do_second_substrate_position(
        &self,
        candidate_a: &mut Genome,
        candidate_b: &mut Genome,
        substrate: &Genome,
        best_substrate: &Genome,
        s1: &Genome,
        s2: &Genome,
        af: f64,
    ) {
        let active_dim = self.params.problem_dimension;
        let ec = self.params.ec;

        let mut rand_vec = vec![0.0; active_dim];
        randomize(&mut rand_vec);

        //  println!("rand_vc : {:?}", rand_vec);

        let sca1: Vec<f64> = rand_vec.iter().map(|x| ec + (1.0 - ec) * x).collect();

        // println!("scA1: {:?}", sca1);

        randomize(&mut rand_vec);

        // println!("rand_vc : {:?}", rand_vec);
        let exa: Vec<f64> = rand_vec
            .iter()
            .map(|x| (ec + (1.0 - ec) * x) * af)
            .collect();

        // println!("exa : {:?}", exa);
        for j in 0..active_dim {
            candidate_a.genes[j] = substrate.genes[j]
                + sca1[j] * (s1.genes[j] - s2.genes[j])
                + (exa[j] * (best_substrate.genes[j] - substrate.genes[j]));
        }

        // Space bound
        for j in 0..active_dim {
            candidate_a.genes[j] = f64::min(candidate_a.genes[j], self.params.upper_bounds[j]);
            candidate_a.genes[j] = f64::max(candidate_a.genes[j], self.params.lower_bounds[j]);
        }
        // ==============================================================================================
        // Scalar random factors  for all dimensions
        let mut rng = rand::thread_rng();
        let intervall = Uniform::from(0.0f64..=1.0);

        let scb1 = ec + (1.0 - ec) * intervall.sample(&mut rng);
        let exb = (ec + (1.0 - ec) * intervall.sample(&mut rng)) * af;
        for j in 0..active_dim {
            candidate_b.genes[j] = substrate.genes[j]
                + scb1 * (s1.genes[j] - s2.genes[j])
                + exb * (best_substrate.genes[j] - substrate.genes[j]);
        }

        // Space bound
        for j in 0..active_dim {
            candidate_b.genes[j] = f64::min(candidate_b.genes[j], self.params.upper_bounds[j]);
            candidate_b.genes[j] = f64::max(candidate_b.genes[j], self.params.lower_bounds[j]);
        }
    }
}

impl<'a, T: Problem> EOA for EAO<'a, T> {
    ///
    /// Call this function to execute GO algorithm.
    ///
    fn run(&mut self) -> OptimizationResult {
        let chronos = Instant::now();

        match self.params.check() {
            Err(error) => OptimizationResult::get_empty(Some(error)),
            Ok(()) => {
                // ------------------------------------------------------------

                #[cfg(feature = "parallel")]
                {
                    let nbr_threads = match self.params.num_threads {
                        None => 1,
                        Some(nthreads) => nthreads,
                    };
                    match rayon::ThreadPoolBuilder::new()
                        .num_threads(nbr_threads)
                        .build_global()
                    {
                        Ok(_) => println!("Thread pool init. for {} threads ... ", nbr_threads),
                        Err(_) => {
                            return OptimizationResult::get_empty(Some(
                                crate::core::OptError::ThreadPoolBuildErr,
                            ))
                        }
                    };
                }
                // ------------------------------------------------------------
                let enzyme_count: usize = self.params.population_size;
                let active_site_dimension: usize = self.params.problem_dimension;

                let max_iter: usize = self.params.max_iterations;
                let mut break_process: bool = false;

                let ub = self.params.upper_bounds;
                let lb = self.params.lower_bounds;

                let mut reaction_rate: Vec<f64> = vec![0.0; enzyme_count];

                let mut convergence_curve = vec![0.0f64; max_iter];

                // random vector of dim length
                let mut rand_vec: Vec<f64> = vec![0.0; active_site_dimension];

                let mut first_substrate_position: Genome =
                    Genome::new(enzyme_count, active_site_dimension);
                let mut second_substrate_position = first_substrate_position.clone();
                let mut updated_position = first_substrate_position.clone();

                let mut _first_evaluation: f64 = 0.0;
                let mut _second_evaluation: f64 = 0.0;
                let mut _updated_fitness: f64 = 0.0;

                // best fitness so-far
                let mut optimal_catalysis: f64 = f64::MAX; // MAX for minimization

                // best solution
                let mut _best_substrate = first_substrate_position.clone();

                let mut s1 = first_substrate_position.clone();
                let mut s2 = first_substrate_position.clone();
                let mut candidate_a = first_substrate_position.clone();
                let mut candidate_b = first_substrate_position.clone();

                let mut _candidate_a_fitness: f64 = 0.0;
                let mut _candidate_b_fitness: f64 = 0.0;
                // =============================================================
                let mut substrate_pool =
                    self.initialize(self.params, InitializationMode::RealUniform);

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
                    for i in 0..enzyme_count {
                        match substrate_pool[i].fitness {
                            None => reaction_rate[i] = f64::MAX,
                            Some(fit) => reaction_rate[i] = fit,
                        };
                    }
                }
                //________________________________________

                // [OptimalCatalysis, idx] = min(ReactionRate);

                let best_idx: usize = reaction_rate
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| a.total_cmp(b))
                    .map(|(i, _)| i)
                    .unwrap_or(0);

                let mut best_substrate = substrate_pool[best_idx].clone();
                // =============================================================

                // MAIN LOOP
                let mut t: usize = 0;
                while t < max_iter {
                    let af = f64::sqrt((t + 1) as f64 / max_iter as f64);

                    for i in 0..enzyme_count {
                        // === 1) Update FirstSubstratePosition

                        // FirstSubstratePosition = (BestSubstrate - SubstratePool(i,:)) +
                        // rand(1, ActiveSiteDimension).* sin(AF * SubstratePool(i,:));
                        randomize(&mut rand_vec);
                        for j in 0..active_site_dimension {
                            first_substrate_position.genes[j] = (best_substrate.genes[j]
                                - substrate_pool[i].genes[j])
                                + (rand_vec[j] * f64::sin(af * substrate_pool[i].genes[j]));
                        }

                        // FirstSubstratePosition = max(min(FirstSubstratePosition, UB),LB);
                        for j in 0..active_site_dimension {
                            first_substrate_position.genes[j] =
                                f64::min(first_substrate_position.genes[j], ub[j]);
                            first_substrate_position.genes[j] =
                                f64::max(first_substrate_position.genes[j], lb[j]);
                        }
                        //  FirstEvaluation = EvaluateCatalysis(FirstSubstratePosition);
                        _first_evaluation = self
                            .problem
                            .objectivefunction(&mut first_substrate_position.genes);

                        // === 2) Pick two random distinct Substrates
                        let (indx1, indx2) = self.pick_two_random_distinct_substrates();
                        // S1 = SubstratePool(Substrates(1), :);
                        // S2 = SubstratePool(Substrates(2), :);
                        copy_vector(&substrate_pool[indx1].genes, &mut s1.genes);
                        copy_vector(&substrate_pool[indx2].genes, &mut s2.genes);
                        // ====== 2.1) vector-valued random factors for each dimension

                        // println!("(s1, s2) = ({}, {})", s1, s2);
                        self.do_second_substrate_position(
                            &mut candidate_a,
                            &mut candidate_b,
                            &substrate_pool[i],
                            &best_substrate,
                            &s1,
                            &s2,
                            af,
                        );

                        _candidate_a_fitness =
                            self.problem.objectivefunction(&mut candidate_a.genes);
                        _candidate_b_fitness =
                            self.problem.objectivefunction(&mut candidate_b.genes);

                        if _candidate_a_fitness < _candidate_b_fitness {
                            copy_solution(
                                &candidate_a,
                                &mut second_substrate_position,
                                active_site_dimension,
                            );
                            _second_evaluation = _candidate_a_fitness;
                        } else {
                            copy_solution(
                                &candidate_b,
                                &mut second_substrate_position,
                                active_site_dimension,
                            );
                            _second_evaluation = _candidate_b_fitness;
                        }

                        //== 3) Compare FirstSubstratePosition vs. SecondSubstratePosition
                        if _first_evaluation < _second_evaluation {
                            copy_solution(
                                &first_substrate_position,
                                &mut updated_position,
                                active_site_dimension,
                            );
                            _updated_fitness = _first_evaluation;
                        } else {
                            copy_solution(
                                &second_substrate_position,
                                &mut updated_position,
                                active_site_dimension,
                            );
                            _updated_fitness = _second_evaluation;
                        };

                        // == 4) Update SubstratePool & Global Best
                        if _updated_fitness < reaction_rate[i] {
                            copy_solution(
                                &updated_position,
                                &mut substrate_pool[i],
                                active_site_dimension,
                            );
                            reaction_rate[i] = _updated_fitness;

                            // Copy the best sol:
                            if _updated_fitness < optimal_catalysis {
                                optimal_catalysis = _updated_fitness;
                                copy_solution(
                                    &updated_position,
                                    &mut best_substrate,
                                    active_site_dimension,
                                );
                            }
                        };
                    }

                    convergence_curve[t] = optimal_catalysis;
                    best_substrate.fitness = Some(optimal_catalysis);

                    self.problem
                        .iteration_increment(t, &best_substrate, &mut break_process);
                    if break_process {
                        t += 1;
                        break;
                    }
                    t += 1;
                }

                best_substrate.fitness =
                    Some(self.problem.objectivefunction(&mut best_substrate.genes));

                let duration = chronos.elapsed();
                let result = OptimizationResult {
                    best_genome: Some(best_substrate),
                    best_fitness: Some(optimal_catalysis),
                    convergence_trend: Some(convergence_curve[0..t].to_vec()),
                    computation_time: Some(duration),
                    err_report: None,
                };
                result
            }
        }
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

    /// Number of threads for parallel execution.
    pub num_threads: Option<usize>,
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
            num_threads: None,
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
            num_threads: None,
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

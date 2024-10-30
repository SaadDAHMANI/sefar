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

    fn clone_population(&mut self, source: &Vec<Genome>, destination: &mut Vec<Genome>) {
        for i in 0..source.len() {
            for j in 0..source[i].genes.len() {
                destination[i].genes[j] = source[i].genes[j];
            }
        }
    }

    fn gained_shared_junior_r1r2r3(
        &self,
        ind_best: &Vec<usize>,
    ) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
        let pop_size = ind_best.len();
        let mut rng = rand::thread_rng();

        // Initialize R1, R2, R3
        let mut r1 = vec![0; pop_size];
        let mut r2 = vec![0; pop_size];

        let interval3 = Uniform::from(1..=pop_size);

        let mut r3: Vec<usize> = (0..pop_size).map(|_| interval3.sample(&mut rng)).collect();

        // R0: Vector from 1 to pop_size
        let r0: Vec<usize> = (1..=pop_size).collect();

        // Fill R1 and R2 according to the position of each element in `ind_best`
        for i in 0..pop_size {
            let ind = ind_best
                .iter()
                .position(|&x| x == i + 1)
                .unwrap_or_default();

            if ind == 0 {
                // Best
                r1[i] = ind_best[1];
                r2[i] = ind_best[2];
            } else if ind == pop_size - 1 {
                // Worst
                r1[i] = ind_best[pop_size - 2];
                r2[i] = ind_best[pop_size - 1];
            } else {
                // Middle
                r1[i] = ind_best[ind - 1];
                r2[i] = ind_best[ind + 1];
            }
        }

        // Generate R3 such that it does not overlap with R1, R2, or R0
        let mut iterations = 0;
        loop {
            let mut conflicts = false;

            for i in 0..pop_size {
                if r3[i] == r1[i] || r3[i] == r2[i] || r3[i] == r0[i] {
                    r3[i] = interval3.sample(&mut rng); //rng.gen_range(1..=pop_size);
                    conflicts = true;
                }
            }

            if !conflicts {
                break;
            }

            iterations += 1;
            if iterations > 1000 {
                //break;
                panic!("Cannot generate R3 without conflicts in 1000 iterations");
            }
        }

        (r1, r2, r3)
    }

    fn gained_shared_senior_r1r2r3(
        &self,
        ind_best: &Vec<usize>,
    ) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
        let pop_size = ind_best.len();

        // Calculate the ranges for R1, R2, and R3
        let r1_size = (pop_size as f64 * 0.1).round() as usize;
        let r2_size = (pop_size as f64 * 0.8).round() as usize;

        let mut rng = rand::thread_rng();

        // R1: First 10% of `ind_best`, then shuffle with random indices
        let r1_slice = &ind_best[0..r1_size];
        let interval_1 = Uniform::from(0..r1_slice.len());

        let mut r1 = Vec::with_capacity(pop_size);
        for _ in 0..pop_size {
            let random_index = interval_1.sample(&mut rng); //rng.gen_range(0..r1_slice.len());
            r1.push(r1_slice[random_index]);
        }

        // R2: Next 80% of `ind_best`, then shuffle with random indices
        let r2_slice = &ind_best[r1_size..r1_size + r2_size];
        let mut r2 = Vec::with_capacity(pop_size);
        let interval_2 = Uniform::from(0..r2_slice.len());

        for _ in 0..pop_size {
            let random_index = interval_2.sample(&mut rng);
            r2.push(r2_slice[random_index]);
        }

        // R3: Last 10% of `ind_best`, then shuffle with random indices
        let r3_slice = &ind_best[r1_size + r2_size..];
        let mut r3 = Vec::with_capacity(pop_size);
        let interval_3 = Uniform::from(0..r3_slice.len());
        for _ in 0..pop_size {
            let random_index = interval_3.sample(&mut rng); //rng.gen_range(0..r3_slice.len());
            r3.push(r3_slice[random_index]);
        }
        (r1, r2, r3)
    }

    fn bound_constraint(&self, vi: &mut Vec<Vec<f64>>, pop: &Vec<Genome>) {
        let np = self.params.population_size; // Population size
        let d = self.params.dimensions; //pop[0].len(); // Dimension
        let lb = self.params.get_lower_bounds();
        let ub = self.params.get_upper_bounds();

        // Check the lower bound
        for i in 0..np {
            for j in 0..d {
                if vi[i][j] < lb[j] {
                    vi[i][j] = (pop[i].genes[j] + lb[j]) / 2.0;
                }
            }
        }
        // Check the upper bound
        for i in 0..np {
            for j in 0..d {
                if vi[i][j] > ub[j] {
                    vi[i][j] = (pop[i].genes[j] + ub[j]) / 2.0;
                }
            }
        }
    }

    fn update_gained_shared_junior_1(
        &self,
        gained_shared_junior: &mut Vec<Vec<f64>>,
        pop: &Vec<Genome>,
        kf: f64,
        ind1: &Vec<bool>,
        rg1: &Vec<usize>,
        rg2: &Vec<usize>,
        rg3: &Vec<usize>,
    ) {
        let problem_size = self.params.dimensions; // pop[0].len();

        for (i, &flag) in ind1.iter().enumerate() {
            if flag {
                //let mut new_row = vec![0.0; problem_size];
                for j in 0..problem_size {
                    gained_shared_junior[i][j] = pop[i].genes[j]
                        + kf * (pop[rg1[i]].genes[j] - pop[rg2[i]].genes[j] + pop[rg3[i]].genes[j]
                            - pop[i].genes[j]);
                }
                //gained_shared_junior[i] = new_row;
            }
        }
    }

    fn update_gained_shared_junior_2(
        &self,
        gained_shared_junior: &mut Vec<Vec<f64>>,
        pop: &Vec<Genome>,
        kf: f64,
        ind1: &Vec<bool>,
        rg1: &Vec<usize>,
        rg2: &Vec<usize>,
        rg3: &Vec<usize>,
    ) {
        let problem_size = self.params.dimensions; // pop[0].len();

        for (i, &flag) in ind1.iter().enumerate() {
            if flag {
                //let mut new_row = vec![0.0; problem_size];
                for j in 0..problem_size {
                    gained_shared_junior[i][j] = pop[i].genes[j]
                        + kf * ((pop[rg1[i]].genes[j] - pop[rg2[i]].genes[j])
                            + (pop[i].genes[j] - pop[rg3[i]].genes[j]));
                }
                //gained_shared_junior[i] = new_row;
            }
        }
    }
}

impl<'a, T: Problem> EOA for GSK<'a, T> {
    fn run(&mut self) -> OptimizationResult {
        let result: OptimizationResult = OptimizationResult::get_none(String::from("n/a"));

        //-------------------------------------------------
        let pop_size: usize = self.params.population_size;
        let max_iter: usize = self.params.max_iterations;
        let problem_size: usize = self.params.dimensions;
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
            println!("fitness[{}] = {}", i, fitness[i]);
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
        let mut g: usize = 0;

        let mut d_gained_shared_junior = vec![0.0f64; pop_size];
        let mut d_gained_shared_senior = vec![0.0f64; pop_size];

        //THE MAIN LOOP
        while nfes < max_nfes {
            g += 1;
            // D_Gained_Shared_Junior=ceil((problem_size)*(1-g/G_Max).^K);
            //   D_Gained_Shared_Senior=problem_size-D_Gained_Shared_Junior;
            for j in 0..pop_size {
                d_gained_shared_junior[j] =
                    problem_size as f64 * (1.0 - g as f64 / g_max as f64).powf(k[j]);
                //println!("d_gained_shared_junior[{}] = {}",j, d_gained_shared_junior[j]);
                d_gained_shared_senior[j] = problem_size as f64 - d_gained_shared_junior[j];
            }

            // clone the old_population to the current one
            self.clone_population(&popold, &mut pop);

            //------------------------------------------------------------
            //Sorte and sorting indexes:
            let mut ind_best: Vec<usize> = (0..fitness.len()).collect();
            ind_best.sort_by(|&a, &b| fitness[a].total_cmp(&fitness[b]));
            //println!("sort indexes are : {:?}", indexes);
            //------------------------------------------------------------

            let (rg1, rg2, rg3) = self.gained_shared_junior_r1r2r3(&ind_best);
            let (r1, r2, r3) = self.gained_shared_senior_r1r2r3(&ind_best);
            //println!("Rg3 : {:?}", rg3);

            // PSEUDO-CODE FOR JUNIOR GAINING SHARING KNOWLEDGE PHASE:
            //Gained_Shared_Junior=zeros(pop_size, problem_size);
            let mut gained_shared_junior = vec![vec![0.0f64; problem_size]; pop_size];
            for i in 0..pop_size {
                if fitness[i] > fitness[rg3[i]] {
                    for j in 0..problem_size {
                        //Gained_Shared_Junior (ind1,:)= pop(ind1,:) + KF*ones(sum(ind1), problem_size) .* (pop(Rg1(ind1),:) - pop(Rg2(ind1),:)+pop(Rg3(ind1), :)-pop(ind1,:)) ;
                        gained_shared_junior[i][j] = pop[i].genes[j]
                            + kf * (pop[rg1[i]].genes[j] - pop[rg2[i]].genes[j]
                                + pop[rg3[i]].genes[j]
                                - pop[i].genes[j]);
                    }
                } else {
                    for j in 0..problem_size {
                        // Gained_Shared_Junior(ind1,:) = pop(ind1,:) + KF*ones(sum(ind1), problem_size) .* (pop(Rg1(ind1),:) - pop(Rg2(ind1),:)+pop(ind1,:)-pop(Rg3(ind1), :)) ;
                        gained_shared_junior[i][j] = pop[i].genes[j]
                            + kf * ((pop[rg1[i]].genes[j] - pop[rg2[i]].genes[j])
                                + (pop[i].genes[j] - pop[rg3[i]].genes[j]));
                    }
                }
            }

            // PSEUDO-CODE FOR SENIOR GAINING SHARING KNOWLEDGE PHASE:
            let mut gained_shared_senior = vec![vec![0.0f64; problem_size]; pop_size];
            for i in 0..pop_size {
                if fitness[i] > fitness[r2[i]] {
                    // Gained_Shared_Senior(ind,:) = pop(ind,:) + KF*ones(sum(ind), problem_size) .* (pop(R1(ind),:) - pop(ind,:) + pop(R2(ind),:) - pop(R3(ind), :)) ;
                    for j in 0..problem_size {
                        gained_shared_senior[i][j] = pop[i].genes[j]
                            + kf * (pop[r1[i]].genes[j] - pop[i].genes[j] + pop[r2[i]].genes[j]
                                - pop[r3[i]].genes[j]);
                    }
                } else {
                    // Gained_Shared_Senior(ind,:) = pop(ind,:) + KF*ones(sum(ind), problem_size) .* (pop(R1(ind),:) - pop(R2(ind),:) + pop(ind,:) - pop(R3(ind), :)) ;
                    for j in 0..problem_size {
                        gained_shared_senior[i][j] = pop[i].genes[j]
                            + kf * (pop[r1[i]].genes[j] - pop[r2[i]].genes[j] + pop[i].genes[j]
                                - pop[r3[i]].genes[j]);
                    }
                }
            }

            // check the lower and the upper bound.
            self.bound_constraint(&mut gained_shared_junior, &pop);
            self.bound_constraint(&mut gained_shared_senior, &pop);

            nfes += 1;
        } // THE MAIN LOOP

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

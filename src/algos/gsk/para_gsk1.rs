use crate::common::*;
use crate::core::eoa::{InitializationMode, EOA};
use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
// use crate::core::OptError;

//use rand::rngs::ThreadRng;
use rand_distr::{Distribution, Uniform};

#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use super::GSKparams;
// use std::fmt::Display;
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
/// To use GSK algorithm:
/// ```rust
/// use sefar::algos::gsk::{GSKparams, GSK};
/// use sefar::benchmarks::functions::Sphere;
/// use sefar::core::optimization_result::OptimizationResult;
/// use sefar::core::eoa::EOA;

/// let settings: GSKparams = GSKparams::default();
/// let mut fo = Sphere {};
/// let mut algo : GSK<Sphere> = GSK::new(&settings, &mut fo);
/// let result: OptimizationResult = algo.run();
/// println!(
///    "Gaining-Sharing Knowledge optimizer (GSK) : F1 (Sphere) test; Result: {:?}",
///    result.to_string());
///```

#[derive(Debug)]
pub struct ParaGSK1<'a, T: Problem> {
    /// The problem to optimize. It must define the Problem trait.
    pub problem: &'a mut T,

    /// Define the parameters of GO algorithm.
    pub params: &'a GSKparams<'a>,
}

impl<'a, T: Problem> ParaGSK1<'a, T> {
    ///
    /// Return a new instance of the Gaining-Sharing Knowledge Optimizer (GSK).
    /// settings: The optimization parameters,
    /// problem: The problem to optimize.
    ///
    /// Example :
    /// ```rust
    /// use sefar::algos::gsk::{GSKparams, GSK};
    /// use sefar::benchmarks::functions::Sphere;
    /// let settings: GSKparams = GSKparams::default();
    /// let mut fo = Sphere {};
    /// let gsk: GSK<Sphere> = GSK::new(&settings, &mut fo);
    /// ```
    pub fn new(settings: &'a GSKparams, problem: &'a mut T) -> Self {
        Self {
            problem,
            params: settings,
        }
    }

    fn evaluate_solutions(&mut self, pop: &mut Vec<Genome>, fitness: &mut [f64]) {
        let n = self.params.get_population_size();
        // Sequential mode

        #[cfg(not(feature = "parallel"))]
        {
            for i in 0..n {
                fitness[i] = self.problem.objectivefunction(&mut pop[i].genes);
                pop[i].fitness = Some(fitness[i]);
                //nfes += 1;
                //println!("fitness[{}] = {}", i, fitness[i]);
            }
        }

        //___________Parallel mode________________
        #[cfg(feature = "parallel")]
        {
            pop.par_iter_mut()
                .for_each(|g| g.fitness = Some(self.problem.objectivefunction(&mut g.genes)));
            for i in 0..n {
                match pop[i].fitness {
                    None => fitness[i] = f64::MAX,
                    Some(fit) => fitness[i] = fit,
                };
            }
        }
        // --------------------------------------
    }

    fn find_indices(&self, x: &Vec<usize>, target: usize) -> usize {
        let y: Vec<usize> = x
            .iter()
            .enumerate()
            .filter_map(|(index, &value)| if value == target { Some(index) } else { None })
            .collect();

        match y.first() {
            Some(index) => *index,
            None => 0usize,
        }
    }

    fn gained_shared_junior_r1r2r3(
        &self,
        ind_best: &Vec<usize>,
        pop_size: usize,
    ) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
        //let pop_size = self.params.population_size;
        let mut rng = rand::thread_rng();

        // Initialize R1, R2, R3
        let mut r1: Vec<usize> = vec![0; pop_size];
        let mut r2: Vec<usize> = vec![0; pop_size];

        let interval3 = Uniform::from(0..pop_size);

        let mut r3: Vec<usize> = (0..pop_size).map(|_| interval3.sample(&mut rng)).collect();

        // R0: Vector from 0 to pop_size-1
        let r0: Vec<usize> = (0..pop_size).collect();

        // Fill R1 and R2 according to the position of each element in `ind_best`
        for i in 0..pop_size {
            let ind = self.find_indices(&ind_best, i);
            if ind == 0 {
                // Best
                r1[i] = ind_best[1];
                r2[i] = ind_best[2];
            } else if ind == pop_size - 1 {
                // Worst
                r1[i] = ind_best[pop_size - 3];
                r2[i] = ind_best[pop_size - 2];
            } else {
                // Middle
                r1[i] = ind_best[ind - 1];
                r2[i] = ind_best[ind + 1];
            }

            //println!("i= {}; ind= {}; R1[i]= {}; R2[i]= {}", i, ind, r1[i], r2[i]);
        }

        //println!("R1 : {:?} \n R2 : {:?}", r1, r2);

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

            if !conflicts || iterations > 1000 {
                break;
            }
            iterations += 1;
        }

        //println!("R1: {:?}, \n R2: {:?}, \n R3: {:?}", r1, r2, r3);

        (r1, r2, r3)
    }

    fn gained_shared_senior_r1r2r3(
        &self,
        ind_best: &Vec<usize>,
        p_ratio: f64,
    ) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
        let pop_size = ind_best.len();
        //let p_ratio = self.params.get_partition_size_p();
        // Calculate the ranges for R1, R2, and R3
        let r1_size = (pop_size as f64 * p_ratio).round() as usize;
        let r2_size = (pop_size as f64 * (1.0 - 2.0 * p_ratio)).round() as usize;

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

    fn bound_constraint(&self, vi: &mut Vec<Vec<f64>>, pop: &Vec<Genome>, np: usize) {
        //let np = self.params.get_population_size(); // Population size
        let d = self.params.problem_dimension; //pop[0].len(); // Dimension
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

    fn generate_d_gained_shared_junior_mask(
        &self,
        d_gained_shared_junior: &Vec<f64>,
        pop_size: usize,
    ) -> Vec<Vec<bool>> {
        // let pop_size: usize = self.params.population_size;
        let problem_size: usize = self.params.problem_dimension;

        // Initialize the mask matrix
        let mut mask = vec![vec![false; problem_size]; pop_size];
        let interval01 = Uniform::from(0.0f64..1.0f64);
        let mut rng = rand::thread_rng();

        for i in 0..pop_size {
            for j in 0..problem_size {
                let random_value: f64 = interval01.sample(&mut rng);
                // Compare random value to (D_Gained_Shared_Junior[i] / problem_size)
                mask[i][j] = random_value <= (d_gained_shared_junior[i] / problem_size as f64);
            }
        }
        mask
    }

    fn generate_d_gained_shared_rand_mask(&self, kr: f64, pop_size: usize) -> Vec<Vec<bool>> {
        //let pop_size: usize = self.params.get_population_size;
        let problem_size: usize = self.params.problem_dimension;
        let interval01 = Uniform::from(0.0f64..1.0f64);
        let mut rng = rand::thread_rng();

        let mut mask = vec![vec![false; problem_size]; pop_size];

        for i in 0..pop_size {
            for j in 0..problem_size {
                let random_value: f64 = interval01.sample(&mut rng);
                mask[i][j] = random_value <= kr;
            }
        }
        mask
    }

    fn and_masks(
        &self,
        mask1: &Vec<Vec<bool>>,
        mask2: &Vec<Vec<bool>>,
        pop_size: usize,
    ) -> Vec<Vec<bool>> {
        //let pop_size: usize = self.params.get_population_size;
        let problem_size: usize = self.params.problem_dimension;

        let mut result_mask = vec![vec![false; problem_size]; pop_size];

        for i in 0..pop_size {
            for j in 0..problem_size {
                result_mask[i][j] = mask1[i][j] && mask2[i][j];
            }
        }
        result_mask
    }

    fn update_gained_shared_junior(
        &self,
        gained_shared_junior: &mut Vec<Vec<f64>>,
        pop: &Vec<Genome>,
        fitness: &Vec<f64>,
        ind_best: &Vec<usize>,
        kf: f64,
        pop_size: usize,
    ) {
        let (rg1, rg2, rg3) = self.gained_shared_junior_r1r2r3(ind_best, pop_size);
        //let pop_size = self.params.population_size;
        let problem_size = self.params.problem_dimension;

        for i in 0..pop_size {
            if fitness[i] > fitness[rg3[i]] {
                for j in 0..problem_size {
                    //Gained_Shared_Junior (ind1,:)= pop(ind1,:) +
                    // KF*ones(sum(ind1), problem_size) .* (pop(Rg1(ind1),:) - pop(Rg2(ind1),:)+
                    // pop(Rg3(ind1), :)-pop(ind1,:)) ;
                    gained_shared_junior[i][j] = pop[i].genes[j]
                        + kf * ((pop[rg1[i]].genes[j] - pop[rg2[i]].genes[j])
                            + (pop[rg3[i]].genes[j] - pop[i].genes[j]));
                }
            } else {
                for j in 0..problem_size {
                    // Gained_Shared_Junior(ind1,:) = pop(ind1,:) +
                    // KF*ones(sum(ind1), problem_size) .* (pop(Rg1(ind1),:)
                    // - pop(Rg2(ind1),:)+pop(ind1,:)-pop(Rg3(ind1), :)) ;
                    gained_shared_junior[i][j] = pop[i].genes[j]
                        + kf * ((pop[rg1[i]].genes[j] - pop[rg2[i]].genes[j])
                            + (pop[i].genes[j] - pop[rg3[i]].genes[j]));
                }
            }
        }
    }

    fn update_gained_shared_senior(
        &self,
        gained_shared_senior: &mut Vec<Vec<f64>>,
        pop: &Vec<Genome>,
        fitness: &Vec<f64>,
        ind_best: &Vec<usize>,
        p: f64,
        kf: f64,
        pop_size: usize,
    ) {
        let (r1, r2, r3) = self.gained_shared_senior_r1r2r3(&ind_best, p);
        //let pop_size = self.params.population_size;
        let problem_size = self.params.problem_dimension;

        for i in 0..pop_size {
            if fitness[i] > fitness[r2[i]] {
                for j in 0..problem_size {
                    // Gained_Shared_Senior(ind,:) = pop(ind,:) +
                    // KF*ones(sum(ind), problem_size) .* (pop(R1(ind),:) - pop(ind,:) +
                    // pop(R2(ind),:) - pop(R3(ind), :)) ;
                    gained_shared_senior[i][j] = pop[i].genes[j]
                        + kf * ((pop[r1[i]].genes[j] - pop[i].genes[j])
                            + (pop[r2[i]].genes[j] - pop[r3[i]].genes[j]));
                }
            } else {
                for j in 0..problem_size {
                    // Gained_Shared_Senior(ind,:) = pop(ind,:) +
                    // KF*ones(sum(ind), problem_size) .* (pop(R1(ind),:) - pop(R2(ind),:) +
                    // pop(ind,:) - pop(R3(ind), :)) ;
                    gained_shared_senior[i][j] = pop[i].genes[j]
                        + kf * ((pop[r1[i]].genes[j] - pop[r2[i]].genes[j])
                            + (pop[i].genes[j] - pop[r3[i]].genes[j]));
                }
            }
        }
    }
}

impl<'a, T: Problem> EOA for ParaGSK1<'a, T> {
    /// Run GSK optimizer:
    /// Example :
    /// ```rust
    /// use sefar::algos::gsk::{GSKparams, GSK};
    /// use sefar::benchmarks::functions::Sphere;
    /// use sefar::core::eoa::EOA;
    ///
    /// let settings: GSKparams = GSKparams::default();
    /// let mut fo = Sphere {};
    /// let mut gsk: GSK<Sphere> = GSK::new(&settings, &mut fo);
    /// let result = gsk.run();
    /// ```
    fn run(&mut self) -> OptimizationResult {
        match Parameters::check(self.params) {
            Err(opt_error) => OptimizationResult::get_empty(Some(opt_error)),
            Ok(()) => {
                match self.params.check() {
                    Err(eror) => OptimizationResult::get_empty(Some(eror)),
                    Ok(()) => {
                        let chronos = Instant::now();

                        //-------------------------------------------------
                        let pop_size: usize = self.params.get_population_size();
                        let max_iter: usize = self.params.max_iterations;
                        let problem_size: usize = self.params.problem_dimension;
                        let mut break_process: bool = false;
                        //let max_nfes: usize = pop_size * (max_iter + 1);
                        //--------------------------------------------------
                        //let mut nfes: usize = 0; // function evaluation counter.
                        let mut bsf_fit_var: f64 = f64::MAX; // the best fitness value.
                        let mut bsf_solution: Genome = Genome::new(pop_size + 1, problem_size); // the best solution
                        let mut fitness: Vec<f64> = vec![0.0f64; pop_size];
                        let mut children_fitness: Vec<f64> = vec![0.0f64; pop_size];
                        let mut run_funcvals: Vec<f64> = vec![0.0f64; max_iter + 1];
                        // ------------------------------------------------------------

                        let mut gained_shared_junior = vec![vec![0.0f64; problem_size]; pop_size];
                        let mut gained_shared_senior = vec![vec![0.0f64; problem_size]; pop_size];
                        // ------------------------------------------------------------

                        let g_max_f64: f64 = max_iter as f64;
                        // Initialize the main population:
                        // Initialize the old population
                        //let mut popold = self.initialize(self.params, InitializationMode::RealUniform);
                        let mut ui = self.initialize(self.params, InitializationMode::RealUniform);

                        // Initialize the current population
                        let mut pop = self.initialize(self.params, InitializationMode::RealUniform); //popold.clone();

                        // let mut objfn_duration: Vec<f64> = vec![0.0; max_iter + 1];

                        // let mut learning_duration: Vec<f64> = vec![0.0; max_iter];
                        // let timer1 = Instant::now();
                        ////  Objective function evaluation:
                        self.evaluate_solutions(&mut pop, &mut fitness);
                        // objfn_duration[0] = timer1.elapsed().as_secs_f64();

                        // Save the best fitness value for convergence trend:
                        for i in 0..pop_size {
                            if fitness[i] < bsf_fit_var {
                                bsf_fit_var = fitness[i];
                                // save the best solution
                                //copy_solution(&pop[i], &mut bsf_solution, problem_size);
                            }
                        }
                        run_funcvals[0] = bsf_fit_var; //save history of convergence.

                        //--------------------------------------------------
                        let p: f64 = self.params.get_partition_size_p();
                        let kf = self.params.kf; //Knowledge Factor.
                        let kr = self.params.kr; //Knowledge Ratio.
                        let k = self.params.k; //Knowledge rate.

                        let mut g: usize = 0;

                        let mut d_gained_shared_junior = vec![0.0f64; pop_size];
                        let mut d_gained_shared_senior = vec![0.0f64; pop_size];

                        let problem_size_f64: f64 = problem_size as f64;

                        //THE MAIN LOOP

                        while g < max_iter {
                            g += 1;
                            // D_Gained_Shared_Junior=ceil((problem_size)*(1-g/G_Max).^K);
                            //   D_Gained_Shared_Senior=problem_size-D_Gained_Shared_Junior;

                            let d_gained_shared_value =
                                problem_size_f64 * ((g_max_f64 - g as f64) / g_max_f64).powf(k);
                            for j in 0..pop_size {
                                d_gained_shared_junior[j] = d_gained_shared_value;
                                //println!("d_gained_shared_junior[{}] = {}",j, d_gained_shared_junior[j]);
                                d_gained_shared_senior[j] =
                                    problem_size_f64 - d_gained_shared_junior[j];
                            }

                            // clone the old_population to the current one
                            // self.clone_population(&popold, &mut pop);
                            /*
                            // Objective function evaluation:
                            for i in 0..pop_size {
                                fitness[i] = self.problem.objectivefunction(&pop[i].genes);
                                pop[i].fitness = Some(fitness[i]);
                                //nfes += 1;
                                //println!("fitness[{}] = {}", i, fitness[i]);
                            }
                            */
                            //------------------------------------------------------------
                            //Sorte and sorting indexes:
                            let mut ind_best: Vec<usize> = (0..fitness.len()).collect();
                            ind_best.sort_by(|&a, &b| fitness[a].total_cmp(&fitness[b]));
                            //println!("fit : {:?} \n sort indexes are : {:?}", fitness, ind_best);
                            //------------------------------------------------------------

                            //let (rg1, rg2, rg3) = self.gained_shared_junior_r1r2r3(&ind_best, pop_size);
                            //println!("Rg3 : {:?}", rg3);
                            //  let (r1, r2, r3) = self.gained_shared_senior_r1r2r3(&ind_best, p);

                            // let learning_timer = Instant::now();

                            // #[cfg(not(feature = "parallel"))]
                            // {
                            //println!("Rg3 : {:?}", rg3);
                            // PSEUDO-CODE FOR JUNIOR GAINING SHARING KNOWLEDGE PHASE:
                            // Gained_Shared_Junior=zeros(pop_size, problem_size);
                            self.update_gained_shared_junior(
                                &mut gained_shared_junior,
                                &pop,
                                &fitness,
                                &ind_best,
                                kf,
                                pop_size,
                            );

                            // PSEUDO-CODE FOR SENIOR GAINING SHARING KNOWLEDGE PHASE:
                            self.update_gained_shared_senior(
                                &mut gained_shared_senior,
                                &pop,
                                &fitness,
                                &ind_best,
                                p,
                                kf,
                                pop_size,
                            );

                            // check the lower and the upper bound.
                            self.bound_constraint(&mut gained_shared_junior, &pop, pop_size);
                            self.bound_constraint(&mut gained_shared_senior, &pop, pop_size);
                            //}

                            //println!("gained_sharied_junior = {:?}", gained_shared_junior);
                            //-------------------------------------------------------------------------------
                            // D_Gained_Shared_Junior_mask=rand(pop_size, problem_size)<=(D_Gained_Shared_Junior(:, ones(1, problem_size))./problem_size);
                            let d_gained_shared_junior_mask = self
                                .generate_d_gained_shared_junior_mask(
                                    &d_gained_shared_junior,
                                    pop_size,
                                );

                            /*println!(
                                "d_gained_shared_junior_mask = {:?}",
                                d_gained_shared_junior_mask
                            );*/

                            //D_Gained_Shared_Senior_mask=~D_Gained_Shared_Junior_mask;
                            let mut d_gained_shared_senior_mask: Vec<Vec<bool>> =
                                vec![vec![false; problem_size]; pop_size];
                            for i in 0..pop_size {
                                for j in 0..problem_size {
                                    d_gained_shared_senior_mask[i][j] =
                                        !d_gained_shared_junior_mask[i][j];
                                }
                            }

                            /*println!(
                                "d_gained_shared_senior_mask = {:?}",
                                d_gained_shared_senior_mask
                            );*/

                            let d_gained_shared_junior_rand_mask =
                                self.generate_d_gained_shared_rand_mask(kr, pop_size);
                            /* println!(
                                "d_gained_shared_junior_rand_mask : {:?}",
                                d_gained_shared_junior_rand_mask
                            ); */

                            let d_gained_shared_junior_mask = self.and_masks(
                                &d_gained_shared_junior_mask,
                                &d_gained_shared_junior_rand_mask,
                                pop_size,
                            );
                            let d_gained_shared_senior_rand_mask =
                                self.generate_d_gained_shared_rand_mask(kr, pop_size);

                            // D_Gained_Shared_Senior_mask=and(D_Gained_Shared_Senior_mask,D_Gained_Shared_Senior_rand_mask);
                            let d_gained_shared_senior_mask = self.and_masks(
                                &d_gained_shared_senior_mask,
                                &d_gained_shared_senior_rand_mask,
                                pop_size,
                            );
                            //ui=pop;

                            for i in 0..pop_size {
                                copy_solution(&pop[i], &mut ui[i], problem_size);
                            }

                            //ui(D_Gained_Shared_Junior_mask) = Gained_Shared_Junior(D_Gained_Shared_Junior_mask);

                            for i in 0..pop_size {
                                for j in 0..problem_size {
                                    if d_gained_shared_junior_mask[i][j] {
                                        ui[i].genes[j] = gained_shared_junior[i][j];
                                    }
                                }
                            }

                            //ui(D_Gained_Shared_Senior_mask) = Gained_Shared_Senior(D_Gained_Shared_Senior_mask);
                            for i in 0..pop_size {
                                for j in 0..problem_size {
                                    if d_gained_shared_senior_mask[i][j] {
                                        ui[i].genes[j] = gained_shared_senior[i][j];
                                    }
                                }
                            }

                            // learning_duration[g - 1] = learning_timer.elapsed().as_secs_f64();

                            //  children_fitness = feval(ui); %
                            /* for i in 0..pop_size {
                                children_fitness[i] = self.problem.objectivefunction(&ui[i].genes);
                                ui[i].fitness = Some(children_fitness[i]);
                                //nfes += 1;
                            }*/
                            // Objective function evaluation for childrens

                            // let timer = Instant::now();

                            self.evaluate_solutions(&mut ui, &mut children_fitness);

                            // objfn_duration[g] = timer.elapsed().as_secs_f64();

                            // SAVE THE BEST SOLUTION:
                            // if children_fitness(i) < bsf_fit_var
                            //    bsf_fit_var = children_fitness(i);
                            //    bsf_solution = ui(i, :);
                            // end
                            for i in 0..pop_size {
                                if children_fitness[i] < bsf_fit_var {
                                    bsf_fit_var = children_fitness[i];
                                    copy_solution(&ui[i], &mut bsf_solution, problem_size);
                                }
                            }
                            /* println!(
                                "iter : {} -- best_fit : {} -- best_sol:{:?}",
                                g, bsf_fit_var, bsf_solution
                            );*/

                            //  #[cfg(feature = "report")]
                            //  println!("Iter : {}, best-fitness : {}", g, bsf_fit_var);

                            // UPDATE THE SEARCH POPULATION:
                            for i in 0..pop_size {
                                if children_fitness[i] < fitness[i] {
                                    //  popold[i] = ui[i].clone();
                                    // copy_solution(&ui[i], &mut popold[i], problem_size);

                                    // COPY BETTER SOULTIONS TO THE SEARCH POPULATION:
                                    copy_solution(&ui[i], &mut pop[i], problem_size);
                                    // COPY THE FITNESS OF THE BETTER SOLUTION TOO:
                                    fitness[i] = children_fitness[i];
                                } /* else {
                                      //popold[i] = pop[i].clone();
                                      copy_solution(&pop[i], &mut popold[i], problem_size);
                                  }*/
                            }

                            // SAVE THE BEST- FITNESS (convergence trend):
                            //run_funcvals = [run_funcvals;bsf_fit_var];
                            run_funcvals[g] = bsf_fit_var;

                            bsf_solution.fitness = Some(bsf_fit_var);

                            self.problem
                                .iteration_increment(g, &bsf_solution, &mut break_process);
                            if break_process {
                                break;
                            }
                        } // THE MAIN LOOP

                        let mut result: OptimizationResult = OptimizationResult::get_empty(None);

                        let duration = chronos.elapsed();
                        /*
                              let objfn_time = objfn_duration.iter().fold(0.0, |sum, t| sum + t);
                              println!("Obj.fun time (total) = {objfn_time} S.");
                              let learn_time = learning_duration.iter().fold(0.0, |sum, t| sum + t);
                              println!("Learning time (total) = {learn_time} S.");
                        */
                        result.best_genome = Some(bsf_solution);
                        result.best_fitness = Some(bsf_fit_var);
                        result.convergence_trend = Some(run_funcvals[0..g + 1].to_vec());
                        result.computation_time = Some(duration);
                        result.err_report = None;
                        result
                    }
                }
            }
        }
    }
}

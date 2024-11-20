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
/// use sefar::algos::gsk::{GSKparams, GSK};
/// use crate::benchmarks::functions::Sphere;
/// let settings: GSKparams = GSKparams::default();
/// let mut fo = Sphere {};
/// let gsk: GSK<Sphere> = GSK::new(&settings, &mut fo);
/// let result: OptimizationResult = algo.run();
/// println!(
///    "Gaining-Sharing Knowledge optimizer (GSK) : F1 (Sphere) test; Result: {:?}",
///    result.to_string());
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
    /// Example :
    /// use sefar::algos::gsk::{GSKparams, GSK};
    /// use crate::benchmarks::functions::Sphere;
    /// let settings: GSKparams = GSKparams::default();
    /// let mut fo = Sphere {};
    /// let gsk: GSK<Sphere> = GSK::new(&settings, &mut fo);
    pub fn new(settings: &'a GSKparams, problem: &'a mut T) -> Self {
        Self {
            problem,
            params: settings,
        }
    }

    fn evaluate_solutions(&mut self, pop: &mut Vec<Genome>, fitness: &mut [f64]) {
        // Sequential mode
        #[cfg(not(feature = "parallel"))]
        {
            for i in 0..self.params.get_population_size() {
                fitness[i] = self.problem.objectivefunction(&pop[i].genes);
                pop[i].fitness = Some(fitness[i]);
                //nfes += 1;
                //println!("fitness[{}] = {}", i, fitness[i]);
            }
        }

        //___________Parallel mode________________
        #[cfg(feature = "parallel")]
        {
            pop.par_iter_mut()
                .for_each(|g| g.fitness = Some(self.problem.objectivefunction(&g.genes)));
            for i in 0..n {
                match x[i].fitness {
                    None => fitness[i] = f64::MAX,
                    Some(fit) => fitness[i] = fit,
                };
            }
        }
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
    ) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
        let pop_size = self.params.population_size;
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

    fn bound_constraint(&self, vi: &mut Vec<Vec<f64>>, pop: &Vec<Genome>) {
        let np = self.params.population_size; // Population size
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
    ) -> Vec<Vec<bool>> {
        let pop_size: usize = self.params.population_size;
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

    fn generate_d_gained_shared_rand_mask(&self, kr: f64) -> Vec<Vec<bool>> {
        let pop_size: usize = self.params.population_size;
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

    fn and_masks(&self, mask1: &Vec<Vec<bool>>, mask2: &Vec<Vec<bool>>) -> Vec<Vec<bool>> {
        let pop_size: usize = self.params.population_size;
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
        rg1: &Vec<usize>,
        rg2: &Vec<usize>,
        rg3: &Vec<usize>,
        kf: f64,
    ) {
        let pop_size = self.params.population_size;
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
        r1: &Vec<usize>,
        r2: &Vec<usize>,
        r3: &Vec<usize>,
        kf: f64,
    ) {
        let pop_size = self.params.population_size;
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

impl<'a, T: Problem> EOA for GSK<'a, T> {
    /// Run GSK optimizer:
    /// Example :
    /// use sefar::algos::gsk::{GSKparams, GSK};
    /// use crate::benchmarks::functions::Sphere;
    /// let settings: GSKparams = GSKparams::default();
    /// let mut fo = Sphere {};
    /// let gsk: GSK<Sphere> = GSK::new(&settings, &mut fo);
    /// let result = gsk.run();
    ///
    fn run(&mut self) -> OptimizationResult {
        match self.params.check() {
            Err(eror) => OptimizationResult::get_none(eror),
            Ok(()) => {
                let chronos = Instant::now();

                //-------------------------------------------------
                let pop_size: usize = self.params.get_population_size();
                let max_iter: usize = self.params.max_iterations;
                let problem_size: usize = self.params.problem_dimension;
                //let max_nfes: usize = pop_size * (max_iter + 1);
                //--------------------------------------------------
                //let mut nfes: usize = 0; // function evaluation counter.
                let mut bsf_fit_var: f64 = f64::MAX; // the best fitness value.
                let mut bsf_solution: Genome = Genome::new(pop_size + 1, problem_size); // the best solution
                let mut fitness: Vec<f64> = vec![0.0f64; pop_size];
                let mut children_fitness: Vec<f64> = vec![0.0f64; pop_size];
                let mut run_funcvals: Vec<f64> = vec![0.0f64; max_iter + 1];
                //--------------------------------------------------

                let g_max_f64: f64 = max_iter as f64;
                // Initialize the main population:
                // Initialize the old population
                //let mut popold = self.initialize(self.params, InitializationMode::RealUniform);
                let mut ui = self.initialize(self.params, InitializationMode::RealUniform);

                // Initialize the current population
                let mut pop = self.initialize(self.params, InitializationMode::RealUniform); //popold.clone();

                // Objective function evaluation:
                self.evaluate_solutions(&mut pop, &mut fitness);

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
                        d_gained_shared_senior[j] = problem_size_f64 - d_gained_shared_junior[j];
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

                    let (rg1, rg2, rg3) = self.gained_shared_junior_r1r2r3(&ind_best);
                    //println!("Rg3 : {:?}", rg3);
                    let (r1, r2, r3) = self.gained_shared_senior_r1r2r3(&ind_best, p);

                    // PSEUDO-CODE FOR JUNIOR GAINING SHARING KNOWLEDGE PHASE:
                    // Gained_Shared_Junior=zeros(pop_size, problem_size);
                    let mut gained_shared_junior = vec![vec![0.0f64; problem_size]; pop_size];
                    self.update_gained_shared_junior(
                        &mut gained_shared_junior,
                        &pop,
                        &fitness,
                        &rg1,
                        &rg2,
                        &rg3,
                        kf,
                    );

                    // PSEUDO-CODE FOR SENIOR GAINING SHARING KNOWLEDGE PHASE:
                    let mut gained_shared_senior = vec![vec![0.0f64; problem_size]; pop_size];
                    self.update_gained_shared_senior(
                        &mut gained_shared_senior,
                        &pop,
                        &fitness,
                        &r1,
                        &r2,
                        &r3,
                        kf,
                    );

                    // check the lower and the upper bound.
                    self.bound_constraint(&mut gained_shared_junior, &pop);
                    self.bound_constraint(&mut gained_shared_senior, &pop);

                    //println!("gained_sharied_junior = {:?}", gained_shared_junior);
                    //-------------------------------------------------------------------------------
                    // D_Gained_Shared_Junior_mask=rand(pop_size, problem_size)<=(D_Gained_Shared_Junior(:, ones(1, problem_size))./problem_size);
                    let d_gained_shared_junior_mask =
                        self.generate_d_gained_shared_junior_mask(&d_gained_shared_junior);

                    /*println!(
                        "d_gained_shared_junior_mask = {:?}",
                        d_gained_shared_junior_mask
                    );*/

                    //D_Gained_Shared_Senior_mask=~D_Gained_Shared_Junior_mask;
                    let mut d_gained_shared_senior_mask: Vec<Vec<bool>> =
                        vec![vec![false; problem_size]; pop_size];
                    for i in 0..pop_size {
                        for j in 0..problem_size {
                            d_gained_shared_senior_mask[i][j] = !d_gained_shared_junior_mask[i][j];
                        }
                    }

                    /*println!(
                        "d_gained_shared_senior_mask = {:?}",
                        d_gained_shared_senior_mask
                    );*/

                    let d_gained_shared_junior_rand_mask =
                        self.generate_d_gained_shared_rand_mask(kr);
                    /* println!(
                        "d_gained_shared_junior_rand_mask : {:?}",
                        d_gained_shared_junior_rand_mask
                    ); */

                    let d_gained_shared_junior_mask = self.and_masks(
                        &d_gained_shared_junior_mask,
                        &d_gained_shared_junior_rand_mask,
                    );
                    let d_gained_shared_senior_rand_mask =
                        self.generate_d_gained_shared_rand_mask(kr);

                    // D_Gained_Shared_Senior_mask=and(D_Gained_Shared_Senior_mask,D_Gained_Shared_Senior_rand_mask);
                    let d_gained_shared_senior_mask = self.and_masks(
                        &d_gained_shared_senior_mask,
                        &d_gained_shared_senior_rand_mask,
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

                    //  children_fitness = feval(ui); %
                    /* for i in 0..pop_size {
                        children_fitness[i] = self.problem.objectivefunction(&ui[i].genes);
                        ui[i].fitness = Some(children_fitness[i]);
                        //nfes += 1;
                    }*/

                    // Objective function evaluation for childrens
                    self.evaluate_solutions(&mut ui, &mut children_fitness);

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
                    // SAVE THE BEST- FITNESS (convergence trend):
                    //run_funcvals = [run_funcvals;bsf_fit_var];
                    run_funcvals[g] = bsf_fit_var;

                    /* println!(
                        "iter : {} -- best_fit : {} -- best_sol:{:?}",
                        g, bsf_fit_var, bsf_solution
                    );*/

                    #[cfg(feature = "report")]
                    println!("Iter : {}, best-fitness : {}", g, bsf_fit_var);

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
                } // THE MAIN LOOP

                let mut result: OptimizationResult =
                    OptimizationResult::get_none(String::from("n/a"));

                let duration = chronos.elapsed();
                result.best_genome = Some(bsf_solution);
                result.best_fitness = Some(bsf_fit_var);
                result.convergence_trend = Some(run_funcvals);
                result.computation_time = Some(duration);
                result.err_report = None;
                result
            }
        }
    }
}

///
/// Define parameters for the Gaining-Sharing Knowledge (GSK) Algorithm.
///
#[derive(Debug, Clone)]
pub struct GSKparams<'a> {
    /// Number of search agents.
    pub population_size: usize,

    /// Dimension of the optimization problem (i.e., the length of the solution).
    pub problem_dimension: usize,

    /// The maximum number of iterations serves as the stopping criterion for the optimization process.
    pub max_iterations: usize,

    /// The lower bounds of the search space.
    pub lower_bounds: &'a [f64],

    /// The upper bounds of the search space.
    pub upper_bounds: &'a [f64],

    /// p is the partition size ratio (p in [0, 1], i.e p varies from 0% to 100%).
    /// Number of best individuals = p*100%;
    /// Number of middle individuals = population_size - 2*p*100%;
    /// Number of worst individuals = p*100%.
    pub partition_size_p: f64,

    /// kf is the knowledge factor parameter (kf > 0). The default value is kf = 0.5.
    pub kf: f64,

    /// kr is the knowledge ratio (kr in [0, 1]). The default value is kr = 0.9.
    pub kr: f64,

    /// k is the knowedge rate (k>0). The default value is k=10.
    pub k: f64,
}

impl<'a> GSKparams<'a> {
    /// Build a new instance of GSKparams, where:
    /// pop_size : population size, i.e., number of search agents;
    /// problem_size : problem dimension, i.e., number of decision variables;
    /// max_iter : maximum number of iterations (stopping criterion);
    /// lb: search space lower bound;
    /// ub: search space upper bound;
    /// partition_size_p : partition size ratio;
    /// kf : knowledge factor;
    /// kr : knowledge ratio;
    /// k : knowledge rate.
    pub fn new(
        pop_size: usize,
        problem_size: usize,
        max_iter: usize,
        lb: &'a [f64],
        ub: &'a [f64],
        p: f64,
        kf: f64,
        kr: f64,
        k: f64,
    ) -> Self {
        Self {
            population_size: pop_size,
            problem_dimension: problem_size,
            max_iterations: max_iter,
            lower_bounds: lb,
            upper_bounds: ub,
            partition_size_p: p,
            kf,
            kr,
            k,
        }
    }

    /// Check the partition size rate 'p'.
    pub fn get_partition_size_p(&self) -> f64 {
        let group1_size: f64 = (self.population_size as f64 * self.partition_size_p).round();
        let group3_size: f64 = (self.population_size as f64 - 2.0 * group1_size).round();
        //println!("group1 : {}, group3: {}", group1_size, group3_size);
        if group1_size < 1.0 || group3_size < 1.0 {
            0.1f64
        } else {
            self.partition_size_p
        }
    }
}

impl<'a> Parameters for GSKparams<'a> {
    fn get_problem_dimension(&self) -> usize {
        self.problem_dimension
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

    fn check(&self) -> Result<(), String> {
        let mut errors: usize = 0;
        let mut msg: String = String::new();

        if self.get_population_size() < 12 {
            msg = String::from("population_size must be greater than or equals 12!; \n");
            errors += 1;
        }

        if self.get_problem_dimension() == 0 {
            msg = format!("{} Search space dimension (i.e., problem dimension or decision variables) must be greater than 0!; \n", msg);
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

        if self.k <= 0.0 {
            msg = format!(
                "{} The knowledge rate 'k' should be greater than 0! [actual value k={:?}]; \n",
                msg, self.k
            );
            errors += 1;
        };

        if self.kf <= 0.0 {
            msg = format!(
                "{} The knowledge factor 'kf' should be greater than 0! [actual value kf={:?}]; \n",
                msg, self.kf
            );
            errors += 1;
        };

        if self.kr < 0.0 || self.kr > 1.0 {
            msg = format!(
                "{} The knowledge ratio 'kr' should be in the range [0, 1]! [actual value kr={:?}]; \n",
                msg, self.kr
            );
            errors += 1;
        };

        if errors > 0 {
            msg = format!("There are [{}] errors : \n {}", errors, msg.trim());
            Err(msg)
        } else {
            Ok(())
        }
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
    ///  GSKparams {
    ///     population_size : 12,
    ///     dimensions : 3,
    ///     max_iterations : 1,
    ///     lower_bounds : &[100.0f64, 100.0, 100.0],
    ///     upper_bounds : &[-100.0f64, -100.0, -100.0],
    ///     p : 0.1,
    ///     kf : 0.5,
    ///     kr : 0.9,
    ///     k : 10.0,
    /// };
    /// ~~~
    ///
    fn default() -> Self {
        GSKparams {
            population_size: 12,
            problem_dimension: 3,
            max_iterations: 1,
            lower_bounds: &[-100.0f64, -100.0, -100.0],
            upper_bounds: &[100.0f64, 100.0, 100.0],
            partition_size_p: 0.1,
            kf: 0.5f64,
            kr: 0.9f64,
            k: 10.0f64,
        }
    }
}

#[cfg(test)]
mod gsk_test {
    use super::*;
    use crate::benchmarks::functions::{Sphere, SumAbsFunction};

    #[test]
    fn gsk_gained_shared_junior_r1r2r3_test_1() {
        let settings: GSKparams = GSKparams::default();
        let mut fo = Sphere {};
        let gsk: GSK<Sphere> = GSK::new(&settings, &mut fo);

        // matlab values : ind_best = [ 6  1  9  8 10  5  7  11  2  4 12  3]
        let ind_best: Vec<usize> = vec![5, 0, 8, 7, 9, 4, 6, 10, 1, 3, 11, 2];

        // matlab values : R1 =[ 6 11 4 2 10  1  5 9  1  8  7  4]
        let ans_r1: Vec<usize> = vec![5, 10, 3, 1, 9, 0, 4, 8, 0, 7, 6, 3];

        // matlab values : R2 = [9 4 12 12 7 9 11 10  8 5 2 3]
        let ans_r2: Vec<usize> = vec![8, 3, 11, 11, 6, 8, 10, 9, 7, 4, 1, 2];

        let (r1, r2, _r3) = gsk.gained_shared_junior_r1r2r3(&ind_best);

        //assert_eq!(gsk.params.population_size, ind_best.len());
        //assert_eq!(r1.len(), ind_best.len());
        assert_eq!(r1, ans_r1);
        assert_eq!(r2, ans_r2);
    }

    #[test]
    fn gsk_gained_shared_junior_r1r2r3_test_2() {
        let mut settings: GSKparams = GSKparams::default();
        settings.population_size = 15;

        let mut fo = Sphere {};
        let gsk: GSK<Sphere> = GSK::new(&settings, &mut fo);

        let ind_best: Vec<usize> = vec![3, 4, 10, 14, 9, 13, 7, 2, 12, 5, 8, 11, 0, 6, 1];

        let ans_r1: Vec<usize> = vec![11, 0, 7, 4, 3, 12, 0, 13, 5, 14, 4, 8, 2, 9, 10];

        let ans_r2: Vec<usize> = vec![6, 6, 12, 10, 10, 8, 1, 2, 11, 13, 14, 0, 5, 7, 9];

        let (r1, r2, _r3) = gsk.gained_shared_junior_r1r2r3(&ind_best);

        assert_eq!(r1, ans_r1);
        assert_eq!(r2, ans_r2);
    }

    #[test]
    fn gsk_gained_shared_senior_r1r2r3_test_1() {
        let settings: GSKparams = GSKparams::default();
        let mut fo = Sphere {};
        let gsk: GSK<Sphere> = GSK::new(&settings, &mut fo);

        let ind_best: Vec<usize> = vec![5, 0, 8, 7, 9, 4, 6, 10, 1, 3, 11, 2];

        let ans_r1: Vec<usize> = vec![5; settings.population_size];

        let ans_r3: Vec<usize> = vec![2; settings.population_size];

        let p = 0.1f64;

        let (r1, r2, r3) = gsk.gained_shared_senior_r1r2r3(&ind_best, p);

        assert_eq!(r1, ans_r1);
        assert_eq!(r3, ans_r3);
        for i in 0..settings.population_size {
            assert_ne!(r1[i], r2[i]);
            assert_ne!(r3[i], r2[i]);
        }
    }

    #[test]
    fn gsk_gained_shared_senior_r1r2r3_test_2() {
        let mut settings: GSKparams = GSKparams::default();
        settings.population_size = 15;

        let mut fo = Sphere {};
        let gsk: GSK<Sphere> = GSK::new(&settings, &mut fo);

        let ind_best: Vec<usize> = vec![3, 4, 10, 14, 9, 13, 7, 2, 12, 5, 8, 11, 0, 6, 1];

        let ans_r3: Vec<usize> = vec![1; settings.population_size];
        let p = 0.1f64;
        let (r1, r2, r3) = gsk.gained_shared_senior_r1r2r3(&ind_best, p);

        let x1: usize = 3;
        let x2: usize = 4;

        assert_eq!(r1.contains(&x1), true);
        assert_eq!(r1.contains(&x2), true);

        for i in 2..ind_best.len() {
            assert_eq!(r1.contains(&ind_best[i]), false);
        }

        assert_eq!(r3, ans_r3);

        for i in 0..settings.population_size {
            assert_ne!(r1[i], r2[i]);
            assert_ne!(r3[i], r2[i]);
        }
    }

    #[test]
    fn gsk_update_gained_shared_junior_test_1() {
        let settings: GSKparams = GSKparams::default();

        let mut fo = SumAbsFunction {};
        let gsk: GSK<SumAbsFunction> = GSK::new(&settings, &mut fo);

        let ind_best: Vec<usize> = vec![8, 2, 0, 11, 9, 4, 3, 6, 5, 7, 1, 10];

        //let ans_r1: Vec<usize> = vec![5, 10, 3, 1, 9, 0, 4, 8, 0, 7, 6, 3];

        //let ans_r2: Vec<usize> = vec![8, 3, 11, 11, 6, 8, 10, 9, 7, 4, 1, 2];

        let (rg1, rg2, _rg3) = gsk.gained_shared_junior_r1r2r3(&ind_best);
        //let rg1 : Vec<usize> = vec![]
        let rg3: Vec<usize> = vec![9, 3, 9, 7, 7, 0, 11, 0, 3, 8, 5, 4];

        let g1: Genome = Genome {
            id: 1,
            genes: vec![-2.7754, 3.7144, -3.2501],
            fitness: None,
        };
        let g2: Genome = Genome {
            id: 2,
            genes: vec![-4.0784, 6.5339, 8.2745],
            fitness: None,
        };
        let g3: Genome = Genome {
            id: 3,
            genes: vec![-5.3083, 2.0748, -1.5022],
            fitness: None,
        };
        let g4: Genome = Genome {
            id: 4,
            genes: vec![5.1778, -1.5360, 8.6996],
            fitness: None,
        };
        let g5: Genome = Genome {
            id: 5,
            genes: vec![-3.2613, -8.8368, -0.6691],
            fitness: None,
        };
        let g6: Genome = Genome {
            id: 6,
            genes: vec![-6.7042, -7.6470, 3.2823],
            fitness: None,
        };

        let g7: Genome = Genome {
            id: 7,
            genes: vec![-4.3057, 4.6932, -7.0230],
            fitness: None,
        };
        let g8: Genome = Genome {
            id: 8,
            genes: vec![-0.8740, 9.1076, -7.8404],
            fitness: None,
        };
        let g9: Genome = Genome {
            id: 9,
            genes: vec![-6.1135, 1.1412, -0.2510],
            fitness: None,
        };

        let g10: Genome = Genome {
            id: 10,
            genes: vec![0.8960, 7.1767, 2.7477],
            fitness: None,
        };
        let g11: Genome = Genome {
            id: 11,
            genes: vec![-4.6326, -9.1553, -5.7893],
            fitness: None,
        };
        let g12: Genome = Genome {
            id: 12,
            genes: vec![-4.1194, 1.3157, 5.3672],
            fitness: None,
        };

        let pop: Vec<Genome> = vec![g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12];
        let mut fitness: Vec<f64> = vec![0.0; settings.population_size];
        for i in 0..settings.population_size {
            fitness[i] = gsk.problem.objectivefunction(&pop[i].genes);
            //pop[i].fitness = Some(fitness[i]);
        }

        assert_eq!(fitness[11], 10.8023);

        let mut gained_shared_junior =
            vec![vec![0.0f64; settings.problem_dimension]; settings.population_size];
        let kf: f64 = 0.5;

        gsk.update_gained_shared_junior(
            &mut gained_shared_junior,
            &pop,
            &fitness,
            &rg1,
            &rg2,
            &rg3,
            kf,
        );

        let ans0: Vec<f64> = vec![-5.205550000000001, 2.3628, -9.6837];
        let ans5 = vec![-6.45565, -4.173500000000001, 0.42479999999999984];
        let ans8 = vec![-13.0256, 1.6600000000000001, -3.8523499999999995];
        let ans11 = vec![-6.38415, 4.6608, 5.386450000000001];

        assert_eq!(gained_shared_junior[0], ans0);
        assert_eq!(gained_shared_junior[5], ans5);
        assert_eq!(gained_shared_junior[8], ans8);
        assert_eq!(gained_shared_junior[11], ans11);
    }

    #[test]
    fn gsk_update_gained_shared_senior_test_1() {
        let settings: GSKparams = GSKparams::default();

        let mut fo = SumAbsFunction {};
        let gsk: GSK<SumAbsFunction> = GSK::new(&settings, &mut fo);

        let ind_best: Vec<usize> = vec![5, 9, 7, 4, 1, 0, 3, 2, 11, 10, 6, 8];

        let (_r1, _r2, _r3) = gsk.gained_shared_senior_r1r2r3(&ind_best, 0.1);

        let r1: Vec<usize> = vec![5; settings.population_size];
        let r3: Vec<usize> = vec![8; settings.population_size];
        let r2: Vec<usize> = vec![2, 7, 2, 4, 6, 10, 7, 6, 9, 0, 3, 1];

        //----check r1 and r3
        assert_eq!(_r1, r1);
        assert_eq!(_r3, r3);
        //-------------------

        let g1: Genome = Genome {
            id: 1,
            genes: vec![9.0244, 3.7681, -3.4864],
            fitness: None,
        };
        let g2: Genome = Genome {
            id: 2,
            genes: vec![9.7028, 1.0832, 4.8482],
            fitness: None,
        };
        let g3: Genome = Genome {
            id: 3,
            genes: vec![5.5498, -4.0074, 7.0592],
            fitness: None,
        };
        let g4: Genome = Genome {
            id: 4,
            genes: vec![9.8544, -2.9158, -3.6107],
            fitness: None,
        };
        let g5: Genome = Genome {
            id: 5,
            genes: vec![-8.0875, -6.3500, 0.2394],
            fitness: None,
        };
        let g6: Genome = Genome {
            id: 6,
            genes: vec![-2.9024, 4.2226, 4.1251],
            fitness: None,
        };

        let g7: Genome = Genome {
            id: 7,
            genes: vec![5.1473, -8.0609, 8.7191],
            fitness: None,
        };
        let g8: Genome = Genome {
            id: 8,
            genes: vec![3.1684, -3.1908, -7.5078],
            fitness: None,
        };
        let g9: Genome = Genome {
            id: 9,
            genes: vec![6.8101, -7.7669, -7.8793],
            fitness: None,
        };

        let g10: Genome = Genome {
            id: 10,
            genes: vec![5.1933, 2.3621, -3.8157],
            fitness: None,
        };
        let g11: Genome = Genome {
            id: 11,
            genes: vec![9.4040, -9.4681, 2.3743],
            fitness: None,
        };
        let g12: Genome = Genome {
            id: 12,
            genes: vec![7.3504, -5.9759, 5.0112],
            fitness: None,
        };

        let pop: Vec<Genome> = vec![g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12];
        let mut fitness: Vec<f64> = vec![0.0; settings.population_size];
        for i in 0..settings.population_size {
            fitness[i] = gsk.problem.objectivefunction(&pop[i].genes);
            //pop[i].fitness = Some(fitness[i]);
        }

        assert_eq!(fitness[11], 18.3375);

        let mut gained_shared_senior =
            vec![vec![0.0f64; settings.problem_dimension]; settings.population_size];
        let kf: f64 = 0.5;

        gsk.update_gained_shared_senior(
            &mut gained_shared_senior,
            &pop,
            &fitness,
            &r1,
            &r2,
            &r3,
            kf,
        );

        let ans: Vec<&[f64]> = vec![
            &[5.905449999999999, 13.6506, -2.7570000000000006],
            &[1.5793499999999998, 4.94095, 4.6724],
            &[0.6935499999999992, 1.9873500000000002, 13.061399999999999],
            &[-3.972800000000001, 1.3618499999999996, 4.316549999999999],
            &[-19.56115, 0.5002000000000004, 2.0017500000000004],
            &[-13.911850000000001, 17.0627, 11.0027],
            &[-0.6984000000000004, 0.3689, 6.607849999999999],
            &[-2.6773000000000007, 5.239000000000001, -9.61905],
            &[1.1454500000000003, 3.292349999999999, 0.15469999999999917],
            &[-1.578500000000001, 7.653849999999999, 2.0218499999999993],
            &[4.77295, -0.19720000000000049, 5.3839999999999995],
            &[3.670349999999999, 3.5484, 10.931899999999999],
        ];
        for i in 0..settings.population_size {
            assert_eq!(gained_shared_senior[i], ans[i]);
        }
    }

    #[test]
    fn gsk_bound_constraint_test_1() {
        let mut junior_set: Vec<Vec<f64>> = vec![
            vec![7.0330e+00, 5.0646e+00, 5.7840e-01],
            vec![-5.1558e+00, -2.4062e+00, 6.3496e+00],
            vec![8.8783e+00, 6.7982e-02, 1.6427e+00],
            vec![1.2844e+01, 2.5021e+00, -1.3194e+01], // of g3
        ];

        let ans_bounded_set: Vec<Vec<f64>> = vec![
            vec![7.033024, 5.064569, 0.578399],
            vec![-5.155761, -2.406182, 6.349578],
            vec![8.878307, 0.067982, 1.642746],
            vec![9.927213, 2.502063, -6.805359],
        ];

        let mut settings: GSKparams = GSKparams::default();
        let lb = vec![-10.0f64; settings.problem_dimension];
        let ub = vec![10.0f64; settings.problem_dimension];
        settings.population_size = 4;
        settings.lower_bounds = &lb;
        settings.upper_bounds = &ub;

        let mut fo = SumAbsFunction {};
        let gsk: GSK<SumAbsFunction> = GSK::new(&settings, &mut fo);

        let g0: Genome = Genome {
            id: 1,
            genes: vec![9.0244, 3.7681, -3.4864],
            fitness: None,
        };
        let g1: Genome = Genome {
            id: 2,
            genes: vec![9.7028, 1.0832, 4.8482],
            fitness: None,
        };
        let g2: Genome = Genome {
            id: 3,
            genes: vec![5.5498, -4.0074, 7.0592],
            fitness: None,
        };
        let g3: Genome = Genome {
            id: 4,
            genes: vec![9.8544, -2.9158, -3.6107],
            fitness: None,
        };

        let pop: Vec<Genome> = vec![g0, g1, g2, g3];

        gsk.bound_constraint(&mut junior_set, &pop);

        for i in 0..settings.population_size {
            for j in 0..settings.problem_dimension {
                assert_eq!(
                    (junior_set[i][j] * 100.0).round(),
                    (ans_bounded_set[i][j] * 100.0).round()
                );
            }
        }
    }

    #[test]
    fn gsk_get_partition_size_p() {
        let mut params: GSKparams = GSKparams::default();
        params.partition_size_p = 1.2;
        assert_eq!(params.get_partition_size_p(), 0.1);

        params.partition_size_p = -0.5;
        assert_eq!(params.get_partition_size_p(), 0.1);
    }
}

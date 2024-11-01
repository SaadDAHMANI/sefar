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

    fn generate_d_gained_shared_junior_mask(
        &self,
        d_gained_shared_junior: &Vec<f64>,
    ) -> Vec<Vec<bool>> {
        let pop_size: usize = self.params.population_size;
        let problem_size: usize = self.params.dimensions;

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

    fn generate_d_gained_shared_junior_rand_mask(&self, kr: f64) -> Vec<Vec<bool>> {
        let pop_size: usize = self.params.population_size;
        let problem_size: usize = self.params.dimensions;
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
        let problem_size: usize = self.params.dimensions;

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
        let problem_size = self.params.dimensions;

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
        let problem_size = self.params.dimensions;

        for i in 0..pop_size {
            if fitness[i] > fitness[r2[i]] {
                for j in 0..problem_size {
                    // Gained_Shared_Senior(ind,:) = pop(ind,:) +
                    // KF*ones(sum(ind), problem_size) .* (pop(R1(ind),:) - pop(ind,:) +
                    // pop(R2(ind),:) - pop(R3(ind), :)) ;
                    gained_shared_senior[i][j] = pop[i].genes[j]
                        + kf * ((pop[r1[i]].genes[j] - pop[r2[i]].genes[j])
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
    fn run(&mut self) -> OptimizationResult {
        let mut result: OptimizationResult = OptimizationResult::get_none(String::from("n/a"));
        let chronos = Instant::now();

        //-------------------------------------------------
        let pop_size: usize = self.params.population_size;
        let max_iter: usize = self.params.max_iterations;
        let problem_size: usize = self.params.dimensions;
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
        let mut popold = self.initialize(self.params, InitializationMode::RealUniform);

        // Initialize the current population
        let mut pop = self.initialize(self.params, InitializationMode::RealUniform); //popold.clone();

        // Objective function evaluation:
        for i in 0..pop_size {
            fitness[i] = self.problem.objectivefunction(&pop[i].genes);
            pop[i].fitness = Some(fitness[i]);
            //nfes += 1;
            //println!("fitness[{}] = {}", i, fitness[i]);
        }
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
        let kf = 0.5; //Knowledge Factor.
        let kr = 0.9; //Knowledge Ratio.
        let k = 10.0;
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
            self.clone_population(&popold, &mut pop);

            //------------------------------------------------------------
            //Sorte and sorting indexes:
            let mut ind_best: Vec<usize> = (0..fitness.len()).collect();
            ind_best.sort_by(|&a, &b| fitness[a].total_cmp(&fitness[b]));
            //println!("fit : {:?} \n sort indexes are : {:?}", fitness, ind_best);
            //------------------------------------------------------------

            let (rg1, rg2, rg3) = self.gained_shared_junior_r1r2r3(&ind_best);
            //println!("Rg3 : {:?}", rg3);
            let (r1, r2, r3) = self.gained_shared_senior_r1r2r3(&ind_best);

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
            //D_Gained_Shared_Senior_mask=~D_Gained_Shared_Junior_mask;
            let mut d_gained_shared_senior_mask: Vec<Vec<bool>> =
                vec![vec![false; problem_size]; pop_size];
            for i in 0..pop_size {
                for j in 0..problem_size {
                    d_gained_shared_senior_mask[i][j] = !d_gained_shared_junior_mask[i][j];
                }
            }

            let d_gained_shared_junior_rand_mask =
                self.generate_d_gained_shared_junior_rand_mask(kr);

            let d_gained_shared_junior_mask = self.and_masks(
                &d_gained_shared_junior_mask,
                &d_gained_shared_junior_rand_mask,
            );
            let d_gained_shared_senior_rand_mask =
                self.generate_d_gained_shared_junior_rand_mask(kr);
            // D_Gained_Shared_Senior_mask=and(D_Gained_Shared_Senior_mask,D_Gained_Shared_Senior_rand_mask);
            let d_gained_shared_senior_mask = self.and_masks(
                &d_gained_shared_senior_mask,
                &d_gained_shared_senior_rand_mask,
            );

            //ui=pop;
            //ui(D_Gained_Shared_Junior_mask) = Gained_Shared_Junior(D_Gained_Shared_Junior_mask);
            let mut ui = pop.clone();
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
            for i in 0..pop_size {
                children_fitness[i] = self.problem.objectivefunction(&ui[i].genes);
                ui[i].fitness = Some(children_fitness[i]);
                //nfes += 1;
            }

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

            println!(
                "iter : {} -- best_fit : {} -- best_sol:{:?}",
                g, bsf_fit_var, bsf_solution
            );

            // UPDATE THE SEARCH POPULATION:
            for i in 0..pop_size {
                if children_fitness[i] < fitness[i] {
                    //  popold[i] = ui[i].clone();
                    copy_solution(&ui[i], &mut popold[i], problem_size);
                } else {
                    //popold[i] = pop[i].clone();
                    copy_solution(&pop[i], &mut popold[i], problem_size);
                }
            }
        } // THE MAIN LOOP

        let duration = chronos.elapsed();
        result.best_genome = Some(bsf_solution);
        result.best_fitness = Some(bsf_fit_var);
        result.convergence_trend = Some(run_funcvals);
        result.computation_time = Some(duration);
        result.err_report = None;
        return result;
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
    ///  GSKparams {
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

#[cfg(test)]
mod gsk_test {
    use crate::benchmarks::functions::Sphere;

    use super::*;

    #[test]
    fn gained_shared_junior_r1r2r3_test_1() {
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
    fn gained_shared_junior_r1r2r3_test_2() {
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
    fn gained_shared_senior_r1r2r3_test_1() {
        let settings: GSKparams = GSKparams::default();
        let mut fo = Sphere {};
        let gsk: GSK<Sphere> = GSK::new(&settings, &mut fo);

        let ind_best: Vec<usize> = vec![5, 0, 8, 7, 9, 4, 6, 10, 1, 3, 11, 2];

        let ans_r1: Vec<usize> = vec![5; settings.population_size];

        let ans_r3: Vec<usize> = vec![2; settings.population_size];

        let (r1, r2, r3) = gsk.gained_shared_senior_r1r2r3(&ind_best);

        assert_eq!(r1, ans_r1);
        assert_eq!(r3, ans_r3);
        for i in 0..settings.population_size {
            assert_ne!(r1[i], r2[i]);
            assert_ne!(r3[i], r2[i]);
        }
    }

    #[test]
    fn gained_shared_senior_r1r2r3_test_2() {
        let mut settings: GSKparams = GSKparams::default();
        settings.population_size = 15;

        let mut fo = Sphere {};
        let gsk: GSK<Sphere> = GSK::new(&settings, &mut fo);

        let ind_best: Vec<usize> = vec![3, 4, 10, 14, 9, 13, 7, 2, 12, 5, 8, 11, 0, 6, 1];

        let ans_r3: Vec<usize> = vec![1; settings.population_size];

        let (r1, r2, r3) = gsk.gained_shared_senior_r1r2r3(&ind_best);

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
    fn update_gained_shared_junior_test_1() {
        let settings: GSKparams = GSKparams::default();

        let mut fo = Sphere {};
        let gsk: GSK<Sphere> = GSK::new(&settings, &mut fo);

        let ind_best: Vec<usize> = vec![8, 2, 0, 11, 9, 4, 3, 6, 5, 7, 1, 10];

        //let ans_r1: Vec<usize> = vec![5, 10, 3, 1, 9, 0, 4, 8, 0, 7, 6, 3];

        //let ans_r2: Vec<usize> = vec![8, 3, 11, 11, 6, 8, 10, 9, 7, 4, 1, 2];

        let (rg1, rg2, _rg3) = gsk.gained_shared_junior_r1r2r3(&ind_best);
        //let rg1 : Vec<usize> = vec![]
        let rg3: Vec<usize> = vec![9, 3, 9, 7, 7, 0, 11, 0, 3, 8, 5, 4];

        let mut g1: Genome = Genome {
            id: 1,
            genes: vec![-2.7754, 3.7144, -3.2501],
            fitness: None,
        };
        let mut g2: Genome = Genome {
            id: 2,
            genes: vec![-4.0784, 6.5339, 8.2745],
            fitness: None,
        };
        let mut g3: Genome = Genome {
            id: 3,
            genes: vec![-5.3083, 2.0748, -1.5022],
            fitness: None,
        };
        let mut g4: Genome = Genome {
            id: 4,
            genes: vec![5.1778, -1.5360, 8.6996],
            fitness: None,
        };
        let mut g5: Genome = Genome {
            id: 5,
            genes: vec![-3.2613, -8.8368, -0.6691],
            fitness: None,
        };
        let mut g6: Genome = Genome {
            id: 6,
            genes: vec![-6.7042, -7.6470, 3.2823],
            fitness: None,
        };

        let mut g7: Genome = Genome {
            id: 7,
            genes: vec![-4.3057, 4.6932, -7.0230],
            fitness: None,
        };
        let mut g8: Genome = Genome {
            id: 8,
            genes: vec![-0.8740, 9.1076, -7.8404],
            fitness: None,
        };
        let mut g9: Genome = Genome {
            id: 9,
            genes: vec![-6.1135, 1.1412, -0.2510],
            fitness: None,
        };

        let mut g10: Genome = Genome {
            id: 10,
            genes: vec![0.8960, 7.1767, 2.7477],
            fitness: None,
        };
        let mut g11: Genome = Genome {
            id: 11,
            genes: vec![-4.6326, -9.1553, -5.7893],
            fitness: None,
        };
        let mut g12: Genome = Genome {
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

        assert_eq!(fitness[11], 10.80);

        let mut gained_shared_junior =
            vec![vec![0.0f64; settings.dimensions]; settings.population_size];
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

        let ans: Vec<f64> = vec![-5.2055e+00, 2.3628e+00, -9.6838e+00];
        assert_eq!(gained_shared_junior[0], ans);
    }
}

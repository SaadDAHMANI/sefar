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

const PI: f64 = std::f64::consts::PI;

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

impl<'a, T: Problem> LshadeSpacma<'a, T> {
    pub fn new(settings: &'a LshadeSpacmaParams, problem: &'a mut T) -> Self {
        Self {
            problem,
            params: settings,
        }
    }
    fn randomize_0to1(&self, randvect: &mut [f64]) {
        let between = Uniform::from(0.0..1.0);
        let mut rng = rand::thread_rng();

        for i in 0..randvect.len() {
            randvect[i] = between.sample(&mut rng);
        }
    }

    fn get_weights(&self, mu: usize) -> Vec<f64> {
        //weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
        let mut weights: Vec<f64> = vec![0.0; mu];
        for i in 1..mu + 1 {
            weights[i - 1] = (mu as f64 + 0.5).ln() - f64::ln(i as f64);
        }

        // weights = weights/sum(weights);
        let sum_weights: f64 = weights.iter().sum();
        for i in 0..mu {
            weights[i] = weights[i] / sum_weights;
        }
        weights
    }

    fn get_mueff(&self, weights: &[f64]) -> f64 {
        //mueff=sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i}
        let sum_weights: f64 = weights.iter().sum();
        let sum_weightsp2: f64 = weights.iter().fold(0.0, |acc, w| acc + w.powi(2));
        sum_weights / sum_weightsp2
    }

    fn strategy_setting_adaptation(
        &self,
        problem_size: usize,
        mueff: f64,
    ) -> (f64, f64, f64, f64, f64) {
        let problem_size_f64 = problem_size as f64;
        //% Strategy parameter setting: Adaptation
        //cc = (4 + mueff/problem_size) / (problem_size + 4 + 2*mueff/problem_size); % time constant for cumulation for C
        let cc = (4.0 + mueff / problem_size_f64)
            / (problem_size_f64 + 4.0 + 2.0 * mueff / problem_size_f64);

        //cs = (mueff+2) / (problem_size+mueff+5);  % t-const for cumulation for sigma control
        let cs = (mueff + 2.0) / (problem_size_f64 + mueff + 5.0);

        //c1 = 2 / ((problem_size+1.3)^2+mueff);    % learning rate for rank-one update of C
        let c1 = 2.0 / ((problem_size_f64 + 1.3).powi(2) + mueff);

        //cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((problem_size+2)^2+mueff));  % and for rank-mu update
        let cmu: f64 = f64::min(
            1.0 - c1,
            2.0 * (mueff - 2.0 + 1.0 / mueff) / ((problem_size_f64 + 2.0).powi(2) + mueff),
        );
        //damps = 1 + 2*max(0, sqrt((mueff-1)/(problem_size+1))-1) + cs; % damping for sigma usually close to 1
        let tmp = f64::sqrt((mueff - 1.0) / (problem_size_f64 + 1.0)) - 1.0;
        let damps = 1.0 + cs + 2.0 * f64::max(0.0, tmp);
        (cc, cs, c1, cmu, damps)
    }

    fn dynamic_strategy_parameters(
        &self,
        problem_size: usize,
    ) -> (
        Vec<f64>,
        Vec<f64>,
        Vec<Vec<usize>>,
        Vec<usize>,
        Vec<Vec<usize>>,
        Vec<Vec<usize>>,
        f64,
    ) {
        // % Initialize dynamic (internal) strategy parameters and constants
        //pc = zeros(problem_size,1);
        //ps = zeros(problem_size,1);   % evolution paths for C and sigma
        //B = eye(problem_size,problem_size);                       % B defines the coordinate system
        //D = ones(problem_size,1);                      % diagonal D defines the scaling
        //C = B * diag(D.^2) * B';            % covariance matrix C
        //invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
        //eigeneval = 0;                      % track update of B and D
        //chiN=problem_size^0.5*(1-1/(4*problem_size)+1/(21*problem_size^2));  % expectation of
        let pc: Vec<f64> = vec![0.0; problem_size];
        let ps: Vec<f64> = vec![0.0; problem_size];
        let mut b: Vec<Vec<usize>> = vec![vec![0; problem_size]; problem_size];
        //let mut c: Vec<Vec<usize>> = vec![vec![0; problem_size]; problem_size];
        for i in 0..problem_size {
            b[i][i] = 1;
            //c[i][i] = 1;
        }

        let d: Vec<usize> = vec![1; problem_size];
        let c = b.clone();
        let invsqrt_c = c.clone();
        let problem_sizef64 = problem_size as f64;
        let chi_n = problem_sizef64.powf(0.5)
            * (1.0 - 1.0 / (4.0 * problem_sizef64) + 1.0 / (21.0 * problem_sizef64.powi(2)));
        (pc, ps, b, d, c, invsqrt_c, chi_n)
    }

    fn clone_population(
        &self,
        source: &Vec<Genome>,
        target: &mut Vec<Genome>,
        problem_size: usize,
    ) {
        for i in 0..source.len() {
            copy_solution(&source[i], &mut target[i], problem_size);
        }
    }

    fn update_memory_params(
        &self,
        pop_size: usize,
        memory_size: usize,
        memory_sf: &[f64],
        memory_cr: &[f64],
    ) -> (Vec<usize>, Vec<f64>, Vec<f64>, Vec<f64>) {
        //mem_rand_index = ceil(memory_size * rand(pop_size, 1));
        let mut rand_vec: Vec<f64> = vec![0.0; pop_size];
        self.randomize_0to1(&mut rand_vec);

        let mut mem_rand_index: Vec<usize> = vec![0; pop_size];
        for i in 0..pop_size {
            mem_rand_index[i] = (rand_vec[i] * memory_size as f64).ceil() as usize - 1;
        }

        println!("mem_rand_index : {:?}", mem_rand_index);

        //mu_sf = memory_sf(mem_rand_index);
        let mu_sf: Vec<f64> = mem_rand_index
            .iter()
            .map(|&index| memory_sf[index])
            .collect();

        //mu_cr = memory_cr(mem_rand_index);
        let mu_cr: Vec<f64> = mem_rand_index
            .iter()
            .map(|&index| memory_cr[index])
            .collect();

        //mem_rand_ratio = rand(pop_size, 1);
        let mut mem_rand_ratio: Vec<f64> = vec![0.0; pop_size];
        self.randomize_0to1(&mut mem_rand_ratio);
        (mem_rand_index, mu_sf, mu_cr, mem_rand_ratio)
    }

    fn generate_crosover_rate(&self, mu_cr: &[f64]) {
        // %% for generating crossover rate
        //cr = normrnd(mu_cr, 0.1);
        //term_pos = find(mu_cr == -1);
        //cr(term_pos) = 0;
        //cr = min(cr, 1);
        //cr = max(cr, 0);
        let mut cr: Vec<f64> = mu_cr.iter().map(|x| normal_rand1(*x, 0.1)).collect();
        for i in 0..cr.len() {
            if cr[i] == -1.0 {
                cr[i] = 0.0;
            }
        }

        let cr: Vec<f64> = cr.iter_mut().map(|x| x.clamp(0.0, 1.0)).collect();
    }

    /// Update the scaling_factor (sf).
    fn update_scaling_factor(
        &self,
        sf: &mut [f64],
        pop_size: usize,
        nfes: usize,
        max_nfes: usize,
        mu_sf: &[f64],
    ) {
        //if(nfes <= max_nfes/2)
        //sf=0.45+.1*rand(pop_size, 1);
        //pos = find(sf <= 0);
        //while ~ isempty(pos)
        //sf(pos)=0.45+0.1*rand(length(pos), 1);
        //pos = find(sf <= 0);
        //end
        //else
        //sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
        //pos = find(sf <= 0);
        //while ~ isempty(pos)
        //sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
        //pos = find(sf <= 0);
        //end
        //end
        //sf = min(sf, 1);

        if nfes < (max_nfes / 2) {
            // sf=0.45+.1*rand(pop_size, 1);
            self.randomize_0to1(sf);
            for sfv in sf.iter_mut() {
                *sfv = (*sfv * 0.1) + 0.45;
            }
        } else {
            //sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
            let interval = Uniform::from(0.0..1.0);
            let mut rng = rand::thread_rng();

            let mut tmp_value: f64 = -1.0;
            for i in 0..pop_size {
                while tmp_value <= 0.0 {
                    tmp_value = mu_sf[i] + 0.1 * (PI * (interval.sample(&mut rng) - 0.5)).tan();
                }
                sf[i] = tmp_value;
            }

            //sf = min(sf, 1);
            for i in 0..pop_size {
                if sf[i] > 1.0 {
                    sf[i] = 1.0
                };
            }
        }
    }

    fn gnr1r2(&self, np1: usize, np2: usize) -> (Vec<usize>, Vec<usize>) {
        // r0 = [1 : pop_size];
        let mut r0: Vec<usize> = vec![0; np1];
        for i in 0..np1 {
            r0[i] = i;
        }
        //
        //NP0 = length(r0);
        //let np0 = r0.len();
        //r1 = floor(rand(1, NP0) * NP1) + 1;
        let mut r1: Vec<usize> = r0.clone();
        let interval = Uniform::from(0..np1);
        let mut rng = rand::thread_rng();
        let mut k: usize = 0;
        for i in 0..np1 {
            while r1[i] == r0[i] && k < 999 {
                r1[i] = interval.sample(&mut rng);
                k += 1;
            }
            k = 0;
        }
        //r2 = floor(rand(1, NP0) * NP2) + 1;
        let mut r2: Vec<usize> = r0.clone();
        let interval2 = Uniform::from(0..np2);
        k = 0;
        for i in 0..np1 {
            while r2[i] == r0[i] || r2[i] == r1[i] {
                r2[i] = interval2.sample(&mut rng);
                k += 1;
                if k > 999 {
                    break;
                }
            }
            k = 0;
        }
        (r1, r2)
    }

    fn choose_from_top_solutions(
        &self,
        pop: &Vec<Genome>,
        pop_size: usize,
        p_best_rate: f64,
        sorted_index: &[usize],
    ) -> Vec<Genome> {
        // pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
        // randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]
        // randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
        // pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions
        //
        let pnpv: f64 = f64::max(p_best_rate * pop_size as f64, 2.0); // choose at least two best solutions
        let pnp: usize = pnpv.round() as usize;

        let mut randindex: Vec<usize> = vec![0; pop_size];
        let interval = Uniform::from(0..pnp);
        let mut rng = rand::thread_rng();
        for i in 0..pop_size {
            randindex[i] = interval.sample(&mut rng);
        }

        let mut pbest: Vec<Genome> = Vec::with_capacity(pop_size);
        for k in 0..pop_size {
            pbest.push(pop[sorted_index[randindex[k]]].clone());
        }
        pbest
    }
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
        // Parameter settings for Hybridization----------------------------
        let first_calss_percentage: f64 = 0.5;
        //-----------------------------------------------------------------
        // Initialize the main population
        let mut pop: Vec<Genome> = self.initialize(self.params, InitializationMode::RealUniform);
        let mut popold: Vec<Genome> = pop.clone();
        //-----------------------------------------------------------------
        let mut fitness: Vec<f64> = vec![0.0; pop_size];
        let mut bsf_fit_var: f64 = f64::MAX;
        let mut bsf_solution: Genome = Genome::new(0, problem_size);
        let mut bsf_index: usize = 0;

        // Fitness function evaluation ------------------------------------
        let mut i: usize = 0;
        for genom in pop.iter_mut() {
            fitness[i] = self.problem.objectivefunction(&genom.genes);
            genom.fitness = Some(fitness[i]);
            i += 1;
            nfes += 1;
        }

        // save the best fitness and the best solution
        for i in 0..pop_size {
            if fitness[i] < bsf_fit_var {
                bsf_fit_var = fitness[i];
                copy_solution(&pop[i], &mut bsf_solution, problem_size);
                bsf_index = i;
            }
        }

        //save the best fitness for convergence trend.
        run_funcvals[0] = bsf_fit_var;

        //----------------------------------------------------------
        //  memory_sf = 0.5 .* ones(memory_size, 1);
        //  memory_cr = 0.5 .* ones(memory_size, 1);
        //  memory_pos = 1;
        let memory_sf = vec![0.5; memory_size];
        let memory_cr = vec![0.5; memory_size];
        let memory_pos = 1; // 0;

        let archive: Archive = Archive {
            np: arc_rate * pop_size as f64,
            pop: vec![0.0; problem_size],
            funvalues: vec![0.0; problem_size],
        };

        let memory_1st_class_percentage: Vec<f64> = vec![first_calss_percentage; memory_size];

        // Initialize CMAES parameters --------------------------------------
        let sigma: f64 = 0.5; // coordinate wise standard deviation (step size)
        let mut xmean: Vec<f64> = vec![0.0; problem_size]; // rand(problem_size, 1); // objective variables initial point
        self.randomize_0to1(&mut xmean);
        let mu: usize = pop_size / 2; //number of parents/points for recombination
        let weights = self.get_weights(mu); //weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
        let mueff: f64 = self.get_mueff(&weights);
        //% Strategy parameter setting: Adaptation----------------------------------------------------------------------
        let (cc, cs, c1, cmu, damps) = self.strategy_setting_adaptation(problem_size, mueff);
        //  % Initialize dynamic (internal) strategy parameters and constants-------------------------------------------
        let (pc, ps, b, d, c, invsqrt_c, chi_n) = self.dynamic_strategy_parameters(problem_size);
        let eigeneval = 0; // track update of B and D
                           //--------------------------------------------------------------------------------------------------------------

        // MAIN LOOP
        let hybridization_flag: usize = 1;

        let mut sf: Vec<f64> = vec![0.0; pop_size];

        while nfes < max_nfes {
            //  pop = popold; the old population becomes the current population
            self.clone_population(&popold, &mut pop, problem_size);
            for i in 0..pop.len() {
                fitness[i] = self.problem.objectivefunction(&pop[i].genes);
                pop[i].fitness = Some(fitness[i]);
                nfes += 1;
            }

            // Sort fitness :
            //Sorte and sorting indexes:
            let mut sorted_index: Vec<usize> = (0..fitness.len()).collect();
            sorted_index.sort_by(|&a, &b| fitness[a].total_cmp(&fitness[b]));

            /* println!(
                "fitness: {:?}, \n sorted_index: {:?}",
                fitness, sorted_index
            );*/
            //mem_rand_index = ceil(memory_size * rand(pop_size, 1));
            //mu_sf = memory_sf(mem_rand_index);
            //mu_cr = memory_cr(mem_rand_index);
            //mem_rand_ratio = rand(pop_size, 1);
            let (mem_rand_index, mu_sf, mu_cr, mem_rand_ratio) =
                self.update_memory_params(pop_size, memory_size, &memory_sf, &memory_cr);
            // Generate crossover rate
            let cr = self.generate_crosover_rate(&mu_cr);

            // Generate scaling factor
            self.update_scaling_factor(&mut sf, pop_size, nfes, max_nfes, &mu_sf);
            //println!("sf = {:?}", sf);

            // For generating Hybridization Class probability
            // Select Class_Select_Index
            let mut class_select_index: Vec<bool> = mem_rand_index
                .iter()
                .map(|&index| {
                    // Adjust MATLAB 1-based index to Rust 0-based index
                    let memory_value = memory_1st_class_percentage[index];
                    memory_value >= mem_rand_ratio[index]
                })
                .collect();

            // Apply Hybridization_flag logic
            if hybridization_flag == 0 {
                // All values will be true (equivalent to MATLAB `or(Class_Select_Index, ~Class_Select_Index)`)
                class_select_index = vec![true; pop_size];
            }

            let mut pop_all = pop.clone(); // Start with pop.
            pop_all.push(Genome::from(0, &archive.pop, f64::MAX)); // Extend with archive.pop
                                                                   //[r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);
            let (r1, r2) = self.gnr1r2(pop_size, pop_all.len());
            // println!("r1 = {:?} \n r2 = {:?}", r1, r2);
            let pbest = self.choose_from_top_solutions(&pop, pop_size, p_best_rate, &sorted_index);
        } //END MAIN LOOP.

        result
    }
}
#[derive(Debug, Clone)]
struct Archive {
    /// The maximum size of the archive.
    pub np: f64,
    /// The solutions stored in te archive.
    pub pop: Vec<f64>,
    /// The function value of the archived solutions.
    pub funvalues: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct LshadeSpacmaParams<'a> {
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

impl<'a> Parameters for LshadeSpacmaParams<'a> {
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
}

impl<'a> Default for LshadeSpacmaParams<'a> {
    /// Return the default values for LshadeSpacmaParams
    /// ```rust
    /// population_size: 10,
    /// problem_dimension: 4,
    /// max_iterations: 2,
    /// lower_bounds: &[-100.0, -100.0, -100.0, -100.0],
    /// upper_bounds: &[100.0, 100.0, 100.0, 100.0],
    /// ```

    fn default() -> Self {
        Self {
            population_size: 10,
            problem_dimension: 4,
            max_iterations: 2,
            lower_bounds: &[-100.0, -100.0, -100.0, -100.0],
            upper_bounds: &[100.0, 100.0, 100.0, 100.0],
        }
    }
}

#[cfg(test)]
mod lshade_spacma_test {
    use crate::benchmarks::functions::Sphere;

    use super::*;

    #[test]
    fn lshade_spacma_get_weights_test1() {
        let settings: LshadeSpacmaParams = LshadeSpacmaParams::default();
        let mut fo: Sphere = Sphere {};
        let algo: LshadeSpacma<Sphere> = LshadeSpacma::new(&settings, &mut fo);
        let mut weights = algo.get_weights(10);
        let ans = vec![
            0.27961, 0.19719, 0.14897, 0.11476, 0.08823, 0.06655, 0.04822, 0.03234, 0.01833, 0.0058,
        ];

        for i in 0..10 {
            weights[i] = (weights[i] * 100000.0).round() / 100000.0;
        }

        let mueff = algo.get_mueff(&weights);

        assert_eq!(weights, ans);
        assert_eq!(mueff, 5.938888539979187);
    }

    #[test]
    fn lshade_spama_strategy_setting_adaptation_test1() {
        let settings: LshadeSpacmaParams = LshadeSpacmaParams::default();
        let mut fo: Sphere = Sphere {};
        let algo: LshadeSpacma<Sphere> = LshadeSpacma::new(&settings, &mut fo);
        let problem_size: usize = 12;
        let mueff: f64 = 5.938804235601242;

        let (cc, cs, c1, cmu, damps) = algo.strategy_setting_adaptation(problem_size, mueff);
        assert_eq!(cc, 0.264564630908058);
        assert_eq!(cs, 0.3460862281251846);
        assert_eq!(c1, 1.093919532188545e-02);
        assert_eq!(cmu, 4.067755393929940e-02);
        assert_eq!(damps, 1.3460862281251846);
    }

    #[test]
    fn lshade_spacma_dynamic_strategy_parameters_test1() {
        let settings: LshadeSpacmaParams = LshadeSpacmaParams::default();
        let mut fo: Sphere = Sphere {};
        let algo: LshadeSpacma<Sphere> = LshadeSpacma::new(&settings, &mut fo);
        let problem_size: usize = 4;

        let (pc, ps, _b, _d, _c, _invsqrt_c, chi_n) =
            algo.dynamic_strategy_parameters(problem_size);
        let ans_pc = vec![0.0f64; 4];
        let ans_ps = vec![0.0f64; 4];
        let ans_chi_n: f64 = 1.880952380952381;
        assert_eq!(pc, ans_pc);
        assert_eq!(ps, ans_ps);
        assert_eq!(chi_n, ans_chi_n);
    }

    #[test]
    fn lshade_spacma_update_memory_params_test1() {
        //mu_sf = memory_sf(mem_rand_index);
        let mem_rand_index = vec![3, 1, 2];
        let memory_sf = vec![10.0, 20.0, 30.0, 40.0];
        let mu_sf: Vec<f64> = mem_rand_index
            .iter()
            .map(|&index| memory_sf[index])
            .collect();

        assert_eq!(mu_sf, vec![40.0, 20.0, 30.0]);
    }

    #[test]
    fn lshade_spacma_crosover_rate_test() {
        let mut cr: Vec<f64> = vec![1.0, -1.0, 0.5, 1.7];
        for i in 0..cr.len() {
            if cr[i] == -1.0 {
                cr[i] = 0.0;
            }
        }

        let cr: Vec<f64> = cr.iter_mut().map(|x| x.clamp(0.0, 1.0)).collect();
        assert_eq!(cr, vec![1.0, 0.0, 0.5, 1.0]);
    }

    #[test]
    fn lshade_spacma_gnr1r2_test() {
        let settings: LshadeSpacmaParams = LshadeSpacmaParams::default();
        let mut fo: Sphere = Sphere {};
        let algo: LshadeSpacma<Sphere> = LshadeSpacma::new(&settings, &mut fo);
        let pop_size = settings.get_population_size();
        let pop_all_size = pop_size + 1;
        let (r1, r2) = algo.gnr1r2(pop_size, pop_all_size);

        for i in 0..pop_size {
            assert_ne!(r1[i], i);
            assert_ne!(r2[i], i);
            assert_ne!(r1[i], r2[i]);
        }
    }
}

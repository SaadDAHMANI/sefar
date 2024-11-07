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

impl<'a, T: Problem> LshadeSpacma<'a, T> {
    pub fn new(settings: &'a LshadeSpacmaParams, problem: &'a mut T) -> Self {
        Self {
            problem,
            params: settings,
        }
    }
    fn randomize_0to1(&self, randvect: &mut Vec<f64>) {
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
        //% Strategy parameter setting: Adaptation---------------------------
        let (cc, cs, c1, cmu, damps) = self.strategy_setting_adaptation(problem_size, mueff);

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
    fn default() -> Self {
        Self {
            population_size: 10,
            problem_dimension: 3,
            max_iterations: 2,
            lower_bounds: &[-100.0, -100.0, -100.0],
            upper_bounds: &[100.0, 100.0, 100.0],
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
}

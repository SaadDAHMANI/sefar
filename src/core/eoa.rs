use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use rand::distributions::{Distribution, Uniform};

///
/// Public trait for Evolutionary Optimization Algorithms
///
pub trait EOA {
    fn initialize<P: Parameters>(&self, params: &P, mode: InitializationMode) -> Vec<Genome> {
        let n: usize = params.get_population_size();
        let dim: usize = params.get_problem_dimension();
        let mut positions: Vec<Genome> = Vec::with_capacity(n);

        let intervall01 = Uniform::from(0.0f64..=1.0f64);
        let mut rng = rand::thread_rng();

        match mode {
            InitializationMode::RealUniform => {
                let lb = &params.get_lower_bounds();
                let ub = &params.get_upper_bounds();

                for i in 0..n {
                    let mut sln = Genome::new(i, dim);

                    for j in 0..dim {
                        //  positions[i][j]= intervall01.sample(&mut rng)*(ub-lb)+lb;
                        sln.genes[j] = intervall01.sample(&mut rng) * (ub[j] - lb[j]) + lb[j];
                    }
                    positions.push(sln);
                }
            }

            InitializationMode::BinaryUnifrom => {
                for i in 0..n {
                    let mut sln = Genome::new(i, dim);
                    for j in 0..dim {
                        if intervall01.sample(&mut rng) < 0.5 {
                            sln.genes[j] = 0.0;
                        } else {
                            sln.genes[j] = 1.0;
                        }
                    }
                    positions.push(sln);
                }
            }
        }

        positions
    }

    ///
    /// Run algorithm until reach stopping criterion and return optiization result
    ///
    fn run(&mut self) -> OptimizationResult {
        OptimizationResult {
            best_genome: None,
            best_fitness: None,
            convergence_trend: None,
            computation_time: None,
            err_report: None,
        }
    }

    ///
    /// Run algorithm until reach stopping criterion and return optiization result
    ///
    fn run_with_params(&mut self, _settings: &impl Parameters) -> OptimizationResult {
        OptimizationResult {
            best_genome: None,
            best_fitness: None,
            convergence_trend: None,
            computation_time: None,
            err_report: None,
        }
    }

    fn randomize(randvect: &mut Vec<f64>) {
        let between = Uniform::from(0.0..=1.0);
        let mut rng = rand::thread_rng();

        for item in randvect.iter_mut() {
            *item = between.sample(&mut rng);
        }
    }
}

pub enum InitializationMode {
    RealUniform,
    BinaryUnifrom,
}

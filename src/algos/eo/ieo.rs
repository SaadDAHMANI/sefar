//
// Implementation of Improved Equilibrium Optimizer (EO)
//

extern crate rand;
use rand::distributions::{Distribution, Uniform};
//use rand::prelude::ThreadRng;
use std::time::Instant;

#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use crate::common::*;
use crate::core::eoa::{InitializationMode, EOA};
use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;

///
/// Equilibrium Optimizer (EO)
/// Reference:
/// "Faramarzi, A., Heidarinejad, M., Stephens, B., & Mirjalili, S. (2020).
/// Equilibrium optimizer: A novel optimization algorithm. Knowledge-Based Systems, 191, 105190."
///
#[derive(Debug)]
pub struct IEO<'a, T: Problem> {
    pub problem: &'a mut T,
    pub params: &'a IEOparams<'a>,
}

impl<'a, T: Problem> IEO<'a, T> {
    pub fn new(settings: &'a IEOparams, problem: &'a mut T) -> Self {
        IEO {
            problem,
            params: settings,
        }
    }
}

impl<'a, T: Problem> EOA for IEO<'a, T> {
    fn run(&mut self) -> OptimizationResult {
        let chronos = Instant::now();

        //check paramaters
        //let params = self.params.clone();

        match self.params.check() {
            Err(error) => OptimizationResult::get_none(error),
            Ok(()) => {
                let dim = self.params.get_problem_dimension();
                let particles_no = self.params.get_population_size();
                let max_iter = self.params.get_max_iterations();

                let ub = self.params.get_upper_bounds();

                let lb = self.params.get_lower_bounds();

                // a1=2;
                // a2=1;
                // GP=0.5;
                let a1: f64 = self.params.a1;
                let a2: f64 = self.params.a2;
                let gp: f64 = self.params.gp;
                let nu: f64 = f64::min(self.params.pool_size_rate, 1.0);
                // Iter=0; V=1;
                let mut iter = 0;
                let v: f64 = 1.0;

                // to store agents fitness values
                let mut ceq1_fit: f64 = f64::MAX;
                //let mut ceq1_index: usize = 0;

                let mut ceq1: Genome = Genome::new(particles_no + 1, dim);

                let mut fitness = vec![0.0f64; particles_no];
                let mut fit_old = vec![0.0f64; particles_no];
                let mut c_old = vec![vec![0.0f64; dim]; particles_no];

                let mut lambda = vec![0.0f64; dim];
                let mut r = vec![0.0f64; dim];
                let mut r1 = vec![0.0f64; dim];
                let mut r2 = vec![0.0f64; dim];
                let mut ceq = vec![0.0f64; dim];
                let mut f = vec![0.0f64; dim];
                let mut _gcp: f64 = 0.0;
                //------------------------------------------
                //let between01 = Uniform::from(0.0..=1.0);
                let mut rng = rand::thread_rng();
                //------------------------------------------

                let mut convergence_curve = vec![0.0f64; max_iter];
                let mut _index: usize = 0;
                let mut _g0: f64 = 0.0;
                let mut _g: f64 = 0.0;

                //C=initialization(Particles_no,dim,ub,lb);

                let mut c: Vec<Genome> =
                    self.initialize(self.params, InitializationMode::RealUniform);

                // the main loop of EO
                while iter < max_iter {
                    // space bound
                    for i in 0..c.len() {
                        for j in 0..dim {
                            if c[i].genes[j] < lb[j] {
                                c[i].genes[j] = lb[j];
                            }

                            if c[i].genes[j] > ub[j] {
                                c[i].genes[j] = ub[j];
                            }
                        }
                    }

                    // compute fitness for search agents
                    // Sequential mode
                    #[cfg(not(feature = "parallel"))]
                    for i in 0..particles_no {
                        fitness[i] = self.problem.objectivefunction(&c[i].genes);
                        c[i].fitness = Some(fitness[i]);
                        //fobj(&c[i]);
                    }

                    // Parallel mode
                    //___________Parallel mode________________
                    #[cfg(feature = "parallel")]
                    {
                        c.par_iter_mut().for_each(|g| {
                            g.fitness = Some(self.problem.objectivefunction(&g.genes))
                        });
                        for i in 0..particles_no {
                            match c[i].fitness {
                                None => fitness[i] = f64::MAX,
                                Some(fit) => fitness[i] = fit,
                            };
                        }
                        //println!("EO: Parallel objective function evaluation was done.");
                    }
                    //________________________________________

                    //-- Memory saving---

                    if iter == 0 {
                        copy_vector(&fitness, &mut fit_old);
                        copy_matrix(&c, &mut c_old);
                    }

                    for i in 0..particles_no {
                        if fit_old[i] < fitness[i] {
                            fitness[i] = fit_old[i];
                            copy_vector2genome(&c_old[i], &mut c[i]);
                        }
                    }

                    copy_matrix(&c, &mut c_old);
                    copy_vector(&fitness, &mut fit_old);

                    //--------------------------- MODIFIED EO -------------------------------

                    // compute the size of the equilibrium pool

                    let jpool: usize =
                        (nu * particles_no as f64 * (1.0 - (iter as f64 / max_iter as f64))).ceil()
                            as usize;

                    // create a matrix to store the pool elements
                    let mut c_pool = vec![vec![0.0f64; dim]; jpool + 1];

                    // interval to choose a random element from the equilibrium pool
                    let interval = Uniform::from(0..jpool);

                    // sort the candidate solutions to choose the equilibrium pool
                    let mut ind: Vec<usize> = (0..particles_no).collect();
                    ind.sort_by(|&a, &b| fitness[a].total_cmp(&fitness[b]));
                    for i in 0..jpool {
                        for j in 0..dim {
                            c_pool[i][j] = c[ind[i]].genes[j];
                        }
                    }
                    // Save the best solution and the best fitness :
                    if ceq1_fit > fitness[ind[0]] {
                        ceq1_fit = fitness[ind[0]];
                        copy_vector2genome(&c[ind[0]].genes, &mut ceq1);
                        ceq1.fitness = Some(ceq1_fit);
                    }

                    //println!("fitness = {:?} \n ind = {:?}", fitness, ind);
                    // compute the average solution
                    let mut sum_value: f64 = 0.0;
                    for j in 0..dim {
                        for i in 0..jpool {
                            sum_value += c_pool[i][j];
                        }
                        c_pool[jpool][j] = sum_value / jpool as f64;
                        sum_value = 0.0;
                    }

                    // comput t using Eq 09
                    let tmpt = (iter / max_iter) as f64;
                    let t: f64 = (1.0 - tmpt).powf(a2 * tmpt);

                    for i in 0..particles_no {
                        randomize(&mut lambda);

                        randomize(&mut r); //  r=rand(1,dim);  r in Eq(11

                        //-------------------------------------------------------
                        // Ceq=C_pool(randi(size(C_pool,1)),:);
                        // random selection of one candidate from the pool
                        _index = interval.sample(&mut rng);
                        copy_vector(&c_pool[_index], &mut ceq);
                        //--------------------------------------------------------
                        // compute F using Eq(11)
                        for j in 0..dim {
                            f[j] = a1
                                * f64::signum(r[j] - 0.5)
                                * (f64::exp(-1.0 * lambda[j] * t) - 1.0);
                        }

                        // r1 and r2 to use them in Eq(15)
                        randomize(&mut r1);
                        randomize(&mut r2);

                        for j in 0..dim {
                            // Eq. 15
                            if r2[j] > gp {
                                _gcp = 0.5 * r1[j];
                            } else {
                                _gcp = 0.0f64;
                            }

                            // Eq. 14
                            _g0 = _gcp * (ceq[j] - lambda[j] * c[i].genes[j]);

                            // Eq 13
                            _g = _g0 * f[j];

                            // Eq. 16
                            c[i].genes[j] = ceq[j]
                                + (c[i].genes[j] - ceq[j]) * f[j]
                                + (_g / (lambda[j] * v)) * (1.0 - f[j]);
                        }
                    }

                    // let duration = chronos.elapsed();
                    // println!("seq--> End computation in : {:?}", duration);

                    convergence_curve[iter] = ceq1_fit;
                    iter += 1;

                    //#[cfg(feature = "report")]
                    println!("Iter : {}, Best-fit : {}", iter, ceq1_fit);
                }

                //return results
                let duration = chronos.elapsed();
                let result = OptimizationResult {
                    best_genome: Some(ceq1),
                    best_fitness: Some(ceq1_fit),
                    convergence_trend: Some(convergence_curve),
                    computation_time: Some(duration),
                    err_report: None,
                };
                return result;
            }
        }
    }
}
/// Define parameters for Improved Equilibrium Optimizer
#[derive(Debug, Clone)]
pub struct IEOparams<'a> {
    /// number of search agents (population size)
    pub population_size: usize,

    /// problem dimension (i.e., number of decision variables)
    pub problem_dimension: usize,

    /// maximum number of iterations
    pub max_iterations: usize,

    /// search space lower bounds
    pub lower_bounds: &'a [f64],

    /// search space upper bounds,
    pub upper_bounds: &'a [f64],

    /// EO parameter
    pub a1: f64,
    /// EO parameter
    pub a2: f64,
    /// EO parameter
    pub gp: f64,
    /// This parameter should be in [0, 1]. It controls the rate of equilibrium pool size reduction.
    pub pool_size_rate: f64,
}

#[allow(dead_code)]
impl<'a> IEOparams<'a> {
    pub fn new(
        p_size: usize,
        dim: usize,
        max_iter: usize,
        lb: &'a [f64],
        ub: &'a [f64],
        a1: f64,
        a2: f64,
        gp: f64,
        pool_size_rate: f64,
    ) -> Result<IEOparams<'a>, String> {
        let params = IEOparams {
            population_size: p_size,
            problem_dimension: dim,
            max_iterations: max_iter,
            lower_bounds: lb,
            upper_bounds: ub,
            a1,
            a2,
            gp,
            pool_size_rate,
        };

        match params.check() {
            Err(error) => Err(error),
            Ok(()) => Ok(params),
        }
    }
}

impl<'a> Parameters for IEOparams<'a> {
    fn get_population_size(&self) -> usize {
        self.population_size
    }

    fn get_problem_dimension(&self) -> usize {
        self.problem_dimension
    }

    fn get_max_iterations(&self) -> usize {
        self.max_iterations
    }

    fn get_lower_bounds(&self) -> Vec<f64> {
        self.lower_bounds.to_vec()
    }

    fn get_upper_bounds(&self) -> Vec<f64> {
        self.upper_bounds.to_vec()
    }
}

impl<'a> Default for IEOparams<'a> {
    ///
    /// Return default values of parameters, as following :
    ///
    /// ~~~
    ///
    ///  use sefar::algos::ieo::*;
    ///
    ///  IEOparams{
    ///     population_size : 10,
    ///     dimensions : 3,
    ///     max_iterations : 100,
    ///     lower_bounds : &[100.0f64, 100.0, 100.0],
    ///     upper_bounds : &[-100.0f64, -100.0, -100.0],
    ///     a1 : 2.0f64,
    ///     a2 : 1.0f64,
    ///     gp : 0.5f64,
    ///     pool_size_rate: 0.2,
    /// };
    /// ~~~
    ///
    fn default() -> Self {
        IEOparams {
            population_size: 10,
            problem_dimension: 3,
            max_iterations: 100,
            lower_bounds: &[-100.0f64, -100.0, -100.0],
            upper_bounds: &[100.0f64, 100.0, 100.0],
            a1: 2.0f64,
            a2: 1.0f64,
            gp: 0.5f64,
            pool_size_rate: 0.2,
        }
    }
}

#[cfg(test)]
mod ieo_params_tests {
    use super::*;

    #[test]
    fn test_ub_slice() {
        let d: usize = 5;
        let n: usize = 10;
        let k: usize = 100;

        let ub = vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let lb = ub.clone();

        let params = IEOparams {
            population_size: n,
            max_iterations: k,
            problem_dimension: d,
            lower_bounds: lb.as_slice(),
            upper_bounds: ub.as_slice(),
            a1: 2.0f64,
            a2: 1.0f64,
            gp: 0.5f64,
            pool_size_rate: 0.0625,
        };

        let sl_ub = vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let slice_ub = sl_ub.as_slice();

        assert_eq!(params.upper_bounds, slice_ub);
    }

    #[test]
    fn test_default_fn() {
        let p = IEOparams::default();
        assert_eq!(p.a1, 2.0f64);
        assert_eq!(p.a2, 1.0f64);
        assert_eq!(p.gp, 0.50f64);
    }

    #[test]
    fn eoparams_unwrap_or_default_test_1() {
        let _ub = vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let _lb = vec![-1.0f64, -2.0, -3.0, -4.0, -5.0];

        let p = IEOparams::new(
            10,
            10,
            100,
            _lb.as_slice(),
            _ub.as_slice(),
            0.5,
            0.5,
            0.5,
            0.2,
        )
        .unwrap_or_default();
        assert_eq!(p.a1, 2.0f64);
        assert_eq!(p.a2, 1.0f64);
        assert_eq!(p.gp, 0.50f64);
    }

    #[test]
    fn eoparams_unwrap_or_default_test_2() {
        let _ub = vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let _lb = vec![-1.0f64, -2.0, -3.0, -4.0, -5.0];

        let p = IEOparams::new(
            10,
            5,
            100,
            _lb.as_slice(),
            _ub.as_slice(),
            0.5,
            0.5,
            0.5,
            0.2,
        )
        .unwrap_or_default();
        assert_eq!(p.a1, 0.50f64);
        assert_eq!(p.a2, 0.50f64);
        assert_eq!(p.gp, 0.50f64);
    }
}

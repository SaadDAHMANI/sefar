//
// Implementation of Binary Equilibrium Optimizer (EO)
//

extern crate rand;
use rand::distributions::{Distribution, Uniform};
use std::time::Instant;

//#[cfg(feature = "parallel")]
//use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use crate::common::*;
use crate::core::eoa::{InitializationMode, EOA};
use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;

///
/// Binary Equilibrium Optimizer (BiEO)
/// Reference:
/// Faramarzi, A., Mirjalili, S., & Heidarinejad, M. (2022).
/// Binary equilibrium optimizer: Theory and application in building optimal control problems.
/// Energy and Buildings, 277, 112503.
/// https://doi.org/10.1016/j.enbuild.2022.112503
///

const A: f64 = 2.0 / std::f64::consts::PI;
const B: f64 = std::f64::consts::PI / 2.0;

#[derive(Debug)]
pub struct BiEO<'a, T: Problem> {
    pub problem: &'a mut T,
    pub params: &'a BiEOparams,
}

impl<'a, T: Problem> BiEO<'a, T> {
    pub fn new(settings: &'a BiEOparams, problem: &'a mut T) -> Self {
        BiEO {
            problem,
            params: settings,
        }
    }
}

impl<'a, T: Problem> EOA for BiEO<'a, T> {
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

                //#[cfg(feature = "binary")]
                let ub: Vec<f64> = vec![1.0; dim];

                //#[cfg(feature = "binary")]
                let lb: Vec<f64> = vec![0.0; dim];

                // a1=2;
                // a2=1;
                // GP=0.5;
                let a1: f64 = self.params.a1;
                let a2: f64 = self.params.a2;
                let gp: f64 = self.params.gp;

                // Initialize variables
                //Ceq1=zeros(1,dim);   Ceq1_fit=inf;
                //Ceq2=zeros(1,dim);   Ceq2_fit=inf;
                //Ceq3=zeros(1,dim);   Ceq3_fit=inf;
                //Ceq4=zeros(1,dim);   Ceq4_fit=inf;

                let mut ceq1 = vec![0.0f64; dim];

                let mut ceq2 = vec![0.0f64; dim];

                let mut ceq3 = vec![0.0f64; dim];

                let mut ceq4 = vec![0.0f64; dim];

                let mut ceq_ave = vec![0.0f64; dim];

                let mut ceq1_fit = f64::MAX;

                let mut ceq2_fit = f64::MAX;

                let mut ceq3_fit = f64::MAX;

                let mut ceq4_fit = f64::MAX;

                let mut ceq1_index: usize = 0;
                //let mut ceq2_index : usize = 0;
                //let mut ceq3_index : usize = 0;
                //let mut ceq4_index : usize = 0;

                // Iter=0; V=1;
                let mut iter = 0;

                // to store agents fitness values
                let mut fitness = vec![0.0f64; particles_no];
                let mut fit_old = vec![0.0f64; particles_no];
                let mut c_old = vec![vec![0.0f64; dim]; particles_no];

                //DeltaC=zeros(Particles_no,dim);
                let mut delta_c = vec![vec![0.0f64; dim]; particles_no];
                // V=zeros(Particles_no,dim);
                let mut v = vec![vec![0.0f64; dim]; particles_no];

                let mut c_pool = vec![vec![0.0f64; dim]; 5];

                let mut lambda = vec![0.0f64; dim];
                let mut r = vec![0.0f64; dim];
                let mut r1 = vec![0.0f64; dim];
                let mut r2 = vec![0.0f64; dim];
                let mut ceq = vec![0.0f64; dim];
                let mut f = vec![0.0f64; dim];
                let mut _gcp: f64 = 0.0;
                //------------------------------------------
                let interval = Uniform::from(0..c_pool.len());
                let interval01 = Uniform::from(0.0..1.0);

                let mut rng = rand::thread_rng();
                //------------------------------------------

                let mut convergence_curve = vec![0.0f64; max_iter];
                let mut _index: usize = 0;
                let mut _g0: f64 = 0.0;
                let mut _g: f64 = 0.0;

                //#[cfg(feature = "binary")]
                let mut c: Vec<Genome> =
                    self.initialize(self.params, InitializationMode::BinaryUnifrom);

                // the main loop of EO
                while iter < max_iter {
                    // compute fitness for search agents
                    // Sequential mode
                    //#[cfg(not(feature = "parallel"))]
                    for i in 0..particles_no {
                        fitness[i] = self.problem.objectivefunction(&c[i].genes);
                        //fobj(&c[i]);
                    }

                    for i in 0..c.len() {
                        // space bound
                        for j in 0..dim {
                            if c[i].genes[j] < lb[j] {
                                c[i].genes[j] = lb[j];
                            }

                            if c[i].genes[j] > ub[j] {
                                c[i].genes[j] = ub[j];
                            }
                        }

                        // check fitness with best
                        if fitness[i] < ceq1_fit {
                            ceq1_index = i;
                            ceq1_fit = fitness[i];
                            copy_vector(&c[i].genes, &mut ceq1);
                        } else if (fitness[i] < ceq2_fit) & (fitness[i] > ceq1_fit) {
                            //ceq2_index = i;
                            ceq2_fit = fitness[i];
                            copy_vector(&c[i].genes, &mut ceq2);
                        } else if (fitness[i] < ceq3_fit)
                            & (fitness[i] > ceq2_fit)
                            & (fitness[i] > ceq1_fit)
                        {
                            //ceq3_index = i;
                            ceq3_fit = fitness[i];
                            copy_vector(&c[i].genes, &mut ceq3);
                        } else if (fitness[i] < ceq4_fit)
                            & (fitness[i] > ceq3_fit)
                            & (fitness[i] > ceq2_fit)
                            & (fitness[i] > ceq1_fit)
                        {
                            //ceq4_index = i;
                            ceq4_fit = fitness[i];
                            copy_vector(&c[i].genes, &mut ceq4);
                        }
                    }

                    //-- Memory saving--------------------------------
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

                    //C_old=C;  fit_old=fitness;
                    copy_matrix(&c, &mut c_old);
                    copy_vector(&fitness, &mut fit_old);
                    // ----------------------------------------------

                    // compute averaged candidate Ceq_ave
                    for i in 0..dim {
                        ceq_ave[i] = ((ceq1[i] + ceq2[i] + ceq3[i] + ceq4[i]) / 4.0).round();
                    }

                    //Equilibrium pool
                    for i in 0..dim {
                        c_pool[0][i] = ceq1[i].round();
                        c_pool[1][i] = ceq2[i].round();
                        c_pool[2][i] = ceq3[i].round();
                        c_pool[3][i] = ceq4[i].round();
                        c_pool[4][i] = ceq_ave[i];
                    }
                    // comput t using Eq 09
                    // t=(1-Iter/Max_iter)^(a2*Iter/Max_iter); % Eq(4)
                    let tmpt = (iter / max_iter) as f64;
                    let t: f64 = (1.0 - tmpt).powf(a2 * tmpt);
                    // let chronos = Instant::now();
                    let mut _alpha: f64 = 0.0;

                    for i in 0..particles_no {
                        randomize(&mut lambda); //  lambda=rand(1,dim);  lambda in Eq(3)
                        randomize(&mut r); //  r=rand(1,dim);

                        //-------------------------------------------------------
                        // Ceq=C_pool(randi(size(C_pool,1)),:);
                        // random selection of one candidate from the pool
                        _index = interval.sample(&mut rng);
                        copy_vector(&c_pool[_index], &mut ceq);
                        //--------------------------------------------------------
                        // compute F using Eq(3)
                        // F=a1*sign(r-0.5).*(exp(-lambda.*t)-1);  %Eq(3)
                        for j in 0..dim {
                            f[j] = a1
                                * f64::signum(r[j] - 0.5)
                                * (f64::exp(-1.0 * lambda[j] * t) - 1.0);
                        }

                        // r1 and r2 to use them in Eq(15)
                        randomize(&mut r1);
                        randomize(&mut r2);

                        //GCP=0.5*rand()*ones(1,dim)*(rand>=GP);  % Eq(7)
                        for j in 0..dim {
                            // Eq. 15
                            if r2[j] > gp {
                                _gcp = 0.5 * r1[j];
                            } else {
                                _gcp = 0.0f64;
                            }

                            // G0=GCP.*(Ceq-lambda.*C(i,:));   % Eq(6)
                            _g0 = _gcp * (ceq[j] - lambda[j] * c[i].genes[j]);

                            // G=G0.*F;   % Eq(5)
                            _g = _g0 * f[j];

                            // Alpha=GP*(rand>GP);   % Eq(13)
                            _alpha = 0.0;
                            if interval01.sample(&mut rng) > gp {
                                _alpha = gp;
                            }
                            // DeltaC(i,:)=((Alpha+C(i,:)-Ceq).*F+(G./lambda).*(1-F)); % Eq(12)
                            delta_c[i][j] = ((_alpha + c[i].genes[j] - ceq[j]) * f[j])
                                + ((_g / lambda[j]) * (1.0 - f[j]));

                            //  V(i,:)=abs((2/pi)*atan((pi/2)*DeltaC(i,:))); % Eq(11)
                            // const A: f64 = 2.0/ std::f64::consts::PI;
                            // const B: f64 = std::f64::consts::PI/2.0;                            //
                            v[i][j] = (A * f64::atan(B * delta_c[i][j])).abs();
                        }

                        //Eq.14
                        if interval01.sample(&mut rng) < 0.5 {
                            for j in 0..dim {
                                if interval01.sample(&mut rng) < v[i][j] {
                                    if c[i].genes[j] == 0.0 {
                                        c[i].genes[j] = 1.0;
                                    } else {
                                        c[i].genes[j] = 0.0;
                                    }
                                }
                            }
                        } else {
                            //Eq.14
                            for j in 0..dim {
                                if interval01.sample(&mut rng) < v[i][j] {
                                    c[i].genes[j] = interval01.sample(&mut rng).round();
                                }
                            }
                        }
                    }

                    convergence_curve[iter] = ceq1_fit;
                    iter += 1;

                    #[cfg(feature = "report")]
                    println!(
                        "Iter : {}, Best-fit : {}, Best-solution : {:?}",
                        iter, ceq1_fit, ceq1
                    );
                }

                //return results
                let duration = chronos.elapsed();
                let result = OptimizationResult {
                    best_genome: Some(Genome::from(ceq1_index, &ceq1, ceq1_fit)),
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
/// Define parameters for Equilibrium Optimizer
#[derive(Debug, Clone)]
pub struct BiEOparams {
    /// number of search agents (population size)
    pub population_size: usize,

    /// problem dimension (i.e., number of decision variables)
    pub problem_dimension: usize,

    /// maximum number of iterations
    pub max_iterations: usize,

    /// BiEO parameter
    pub a1: f64,
    /// EO parameter
    pub a2: f64,
    /// EO parameter
    pub gp: f64,
}

#[allow(dead_code)]
impl BiEOparams {
    pub fn new(
        p_size: usize,
        dim: usize,
        max_iter: usize,
        a1: f64,
        a2: f64,
        gp: f64,
    ) -> Result<BiEOparams, String> {
        let params = BiEOparams {
            population_size: p_size,
            problem_dimension: dim,
            max_iterations: max_iter,
            a1,
            a2,
            gp,
        };

        match params.check() {
            Err(error) => Err(error),
            Ok(()) => Ok(params),
        }
    }
}

impl Parameters for BiEOparams {
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
        let lb: Vec<f64> = vec![0.0; self.get_problem_dimension()];
        lb
    }

    fn get_upper_bounds(&self) -> Vec<f64> {
        let ub: Vec<f64> = vec![1.0; self.get_problem_dimension()];
        ub
    }
}

impl Default for BiEOparams {
    ///
    /// Return default values of parameters, as following :
    ///
    /// ~~~
    ///
    ///  use sefar::algos::eo::*;
    ///
    ///  EOparams{
    ///     population_size : 10,
    ///     dimensions : 3,
    ///     max_iterations : 100,
    ///     a1 : 2.0f64,
    ///     a2 : 1.0f64,
    ///     gp : 0.5f64,
    /// };
    /// ~~~
    ///
    fn default() -> Self {
        BiEOparams {
            population_size: 10,
            problem_dimension: 3,
            max_iterations: 100,
            a1: 2.0f64,
            a2: 1.0f64,
            gp: 0.5f64,
        }
    }
}

#[cfg(test)]
mod eo_params_tests {
    use super::*;

    #[test]
    fn test_default_fn() {
        let p = BiEOparams::default();
        assert_eq!(p.a1, 2.0f64);
        assert_eq!(p.a2, 1.0f64);
        assert_eq!(p.gp, 0.50f64);
    }
}

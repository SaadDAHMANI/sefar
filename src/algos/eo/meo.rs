//
// Implementation of Modified Equilibrium Optimizer (m-EO)
//

extern crate rand;
use rand::distributions::{Distribution, Uniform};
//use rand::prelude::ThreadRng;
use std::time::Instant;

use crate::algos::eo::EOparams;
use crate::common::*;
use crate::core::eoa::{InitializationMode, EOA};
use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;

///
/// Sequential Modified Equilibrium Optimizer (MEO)
/// Reference:
/// "Gupta, S., Deep, K., & Mirjalili, S. (2020).
/// An efficient equilibrium optimizer with mutation strategy for numerical optimization.
/// Applied Soft Computing, 96, 106542."
///
#[derive(Debug)]
pub struct MEO<'a, T: Problem> {
    pub problem: &'a mut T,
    pub params: &'a EOparams<'a>,
    pub optimization_result: OptimizationResult,
}

impl<'a, T: Problem> MEO<'a, T> {
    pub fn new(settings: &'a EOparams, problem: &'a mut T) -> Self {
        let result = OptimizationResult {
            best_genome: None,
            best_fitness: None,
            convergence_trend: None,
            computation_time: None,
            err_report: None,
        };

        MEO {
            problem,
            params: settings,
            optimization_result: result,
        }
    }
}

impl<'a, T: Problem> EOA for MEO<'a, T> {
    fn run(&mut self) -> OptimizationResult {
        let chronos = Instant::now();

        //check paramaters
        //let params = self.params.clone();

        match self.params.check() {
            Err(error) => OptimizationResult::get_empty(Some(error)),
            Ok(()) => {
                let dim = self.params.get_problem_dimension(); //self.params.get_dimensions();
                let particles_no = self.params.get_population_size(); //self.params.get_population_size();
                let lb = self.params.get_lower_bounds(); //self.params.get_lower_bounds();
                let ub = self.params.get_upper_bounds();
                let max_iter = self.params.get_max_iterations();

                let mut break_process: bool = false;
                //
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

                // Iter=0; V=1;
                let mut iter = 0;
                let v: f64 = 1.0;

                // to store agents fitness values
                let mut fitness = vec![0.0f64; particles_no];
                let mut fit_old = vec![0.0f64; particles_no];
                let mut c_old = vec![vec![0.0f64; dim]; particles_no];
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
                //let between01 = Uniform::from(0.0..=1.0);
                let mut rng = rand::thread_rng();
                //------------------------------------------

                let mut convergence_curve = vec![0.0f64; max_iter];
                let mut _index: usize = 0;
                let mut _g0: f64 = 0.0;
                let mut _g: f64 = 0.0;

                //let chronos = Instant::now();

                // Step 1: initialize the population randomly within the solution space
                let mut c = self.initialize(self.params, InitializationMode::RealUniform);

                // Step 2 : Evaluate the fitness value of each candidate soluion
                for genom in c.iter_mut() {
                    genom.fitness = Some(self.problem.objectivefunction(&mut genom.genes));
                }

                // Step 3 : Select the 4 best solutions
                // 3.1 Sorting :
                c.sort_by(Genome::cmp_genome);

                // 3.2 update indexes
                for i in 0..c.len() {
                    c[i].id = i;
                }

                // check
                /*
                                for g in c.iter() {
                                    println!("id: {}, fit: {:?}", g.id, g.fitness);
                                }
                */
                // the main loop of EO
                while iter < max_iter {
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

                        // compute fitness for agents

                        fitness[i] = self.problem.objectivefunction(&mut c[i].genes); //fobj(&c[i]);

                        // check fitness with best
                        if fitness[i] < ceq1_fit {
                            ceq1_index = i;
                            ceq1_fit = fitness[i];
                            //copy_vector(&c[i].genes, &mut ceq1);
                            ceq1[..dim].clone_from_slice(&c[i].genes[..dim]);
                        } else if (fitness[i] < ceq2_fit) & (fitness[i] > ceq1_fit) {
                            //ceq2_index = i;
                            ceq2_fit = fitness[i];
                            //copy_vector(&c[i].genes, &mut ceq2);
                            ceq2[..dim].clone_from_slice(&c[i].genes[..dim]);
                        } else if (fitness[i] < ceq3_fit)
                            & (fitness[i] > ceq2_fit)
                            & (fitness[i] > ceq1_fit)
                        {
                            //ceq3_index = i;
                            ceq3_fit = fitness[i];
                            //copy_vector(&c[i].genes, &mut ceq3);
                            ceq3[..dim].clone_from_slice(&c[i].genes[..dim]);
                        } else if (fitness[i] < ceq4_fit)
                            & (fitness[i] > ceq3_fit)
                            & (fitness[i] > ceq2_fit)
                            & (fitness[i] > ceq1_fit)
                        {
                            //ceq4_index = i;
                            ceq4_fit = fitness[i];
                            //copy_vector(&c[i].genes, &mut ceq4);
                            ceq4[..dim].clone_from_slice(&c[i].genes[..dim]);
                        }
                    }

                    // copy the best 4 genomes
                    //copy_vector(&c[ceq1_index].genes, &mut ceq1);
                    //copy_vector(&c[ceq2_index].genes, &mut ceq2);
                    //copy_vector(&c[ceq3_index].genes, &mut ceq3);
                    //copy_vector(&c[ceq4_index].genes, &mut ceq4);

                    //ceq1_fit = fitness[ceq1_index];
                    //ceq2_fit = fitness[ceq2_index];
                    //ceq3_fit = fitness[ceq3_index];
                    //ceq4_fit = fitness[ceq4_index];

                    //-- Memory saving---

                    if iter == 0 {
                        //copy_vector(&fitness, &mut fit_old);
                        fit_old[..particles_no].clone_from_slice(&fitness[..particles_no]);
                        copy_matrix(&c, &mut c_old);
                    }

                    for i in 0..particles_no {
                        if fit_old[i] < fitness[i] {
                            fitness[i] = fit_old[i];
                            //copy_vector2genome(&c_old[i], &mut c[i]);
                            c[i].genes[..dim].clone_from_slice(&c_old[i][..dim]);
                        }
                    }

                    copy_matrix(&c, &mut c_old);
                    //copy_vector(&fitness, &mut fit_old);
                    fit_old[..particles_no].clone_from_slice(&fitness[..particles_no]);
                    // compute averaged candidate Ceq_ave
                    for i in 0..dim {
                        ceq_ave[i] = (ceq1[i] + ceq2[i] + ceq3[i] + ceq4[i]) / 4.0;
                    }

                    //Equilibrium pool
                    c_pool[0][..dim].clone_from_slice(&ceq1[..dim]);
                    c_pool[1][..dim].clone_from_slice(&ceq2[..dim]);
                    c_pool[2][..dim].clone_from_slice(&ceq3[..dim]);
                    c_pool[3][..dim].clone_from_slice(&ceq4[..dim]);
                    c_pool[4][..dim].clone_from_slice(&ceq_ave[..dim]);

                    // comput t using Eq 09
                    let tmpt = (iter / max_iter) as f64;
                    let t: f64 = (1.0 - tmpt).powf(a2 * tmpt);

                    // let chronos = Instant::now();

                    for i in 0..particles_no {
                        MEO::<'a, T>::randomize(&mut lambda); //  lambda=rand(1,dim);  lambda in Eq(11)
                        MEO::<'a, T>::randomize(&mut r); //  r=rand(1,dim);  r in Eq(11

                        //-------------------------------------------------------
                        // Ceq=C_pool(randi(size(C_pool,1)),:);
                        // random selection of one candidate from the pool
                        _index = interval.sample(&mut rng);

                        //copy_vector(&c_pool[_index], &mut ceq);
                        ceq[..dim].clone_from_slice(&c_pool[_index][..dim]);
                        //--------------------------------------------------------
                        // compute F using Eq(11)
                        for j in 0..dim {
                            f[j] = a1
                                * f64::signum(r[j] - 0.5)
                                * (f64::exp(-1.0 * lambda[j] * t) - 1.0);
                        }

                        // r1 and r2 to use them in Eq(15)
                        MEO::<'a, T>::randomize(&mut r1);
                        MEO::<'a, T>::randomize(&mut r2);

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

                    self.problem.iteration_increment(
                        iter,
                        &Genome::from(ceq1_index, &ceq1, ceq1_fit),
                        &mut break_process,
                    );
                    if break_process {
                        break;
                    }
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
                // copy result to EO struct
                self.optimization_result = result.clone();
                result
            }
        }
    }
}

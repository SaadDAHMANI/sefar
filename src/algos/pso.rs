//
// Implementation of Particle Swarm Optimization algorithm (PSO)
//

extern crate rand;
use std::fmt::Display;
//use rand::distributions::{Distribution, Uniform};
//use rand::prelude::ThreadRng;
use std::time::Instant;

use crate::core::eoa::{InitializationMode, EOA};
use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
//use crate::common::*;

pub struct PSO<'a, T: Problem> {
    pub problem: &'a mut T,
    pub params: &'a PSOparams<'a>,
}

impl<'a, T: Problem> PSO<'a, T> {
    pub fn new(settings: &'a PSOparams, problem: &'a mut T) -> Self {
        PSO {
            problem,
            params: settings,
        }
    }
}

impl<'a, T: Problem> EOA for PSO<'a, T> {
    fn run(&mut self) -> OptimizationResult {
        // start time computation
        let chronos = Instant::now();

        //check paramaters
        //let params = self.clone();

        match self.params.check() {
            Err(error) => OptimizationResult::get_none(error),
            Ok(()) => {
                let dim = self.params.problem_dimension;
                let ub = self.params.upper_bounds;
                let lb = self.params.lower_bounds;
                let max_iter = self.params.max_iterations;
                let nop = self.params.population_size;

                // Define the PSO's paramters
                let c1: f64 = self.params.c1;
                let c2: f64 = self.params.c2;

                let w_max: f64 = 0.9;
                let w_min: f64 = 0.2;

                let mut v_max = Vec::new();
                let mut v_min = Vec::new();

                for i in 0..dim {
                    v_max.push((ub[i] - lb[i]) * 0.2f64);
                    v_min.push(-1.0 * v_max[i]);
                }

                let mut cgcurve = vec![0.0f64; max_iter];

                // Velocities initialization
                let mut v = vec![vec![0.0f64; dim]; nop];

                //let mut currentx = Solution::new(nop+1, dim);

                let mut gbest_x = Genome::new(dim + 1, dim);
                //let mut gbest_0 = Vec::new();

                let mut rand1 = vec![0.0f64; dim];
                let mut rand2 = vec![0.0f64; dim];

                // PSO algorithm
                // Particles initialization
                let mut particles = self.initialize(self.params, InitializationMode::RealUniform);

                //initialize pbest_x population with (fitness = f64::MAX)
                let mut pbest_x = self.initialize(self.params, InitializationMode::RealUniform);

                // Main PSO loop
                for t in 0..max_iter {
                    //let mut gbest_index : usize = 0;

                    for k in 0..nop {
                        //Objective function computation
                        // Evaluate search agent using objective function
                        particles[k].fitness =
                            Some(self.problem.objectivefunction(&mut particles[k].genes));

                        //Update the pbest
                        if particles[k].fitness < pbest_x[k].fitness {
                            //pbest_x[k] = particles[k].clone();
                            //copy_genome(&particles[k], &mut pbest_x[k]);
                            pbest_x[k].genes[..dim].clone_from_slice(&particles[k].genes[..dim]);
                        }

                        //Update the gbest
                        if particles[k].fitness < gbest_x.fitness {
                            // gbest_index = k;
                            gbest_x = particles[k].clone();
                            //copy_genome(&particles[k], &mut gbest_x);
                        }
                    }

                    //Update the x and v
                    let tf64 = t as f64;
                    let max_iterf64 = max_iter as f64;
                    let w = w_max - ((tf64 * (w_max - w_min)) / max_iterf64);

                    for k in 0..nop {
                        PSO::<'a, T>::randomize(&mut rand1);
                        PSO::<'a, T>::randomize(&mut rand2);

                        for j in 0..dim {
                            v[k][j] = (w * v[k][j])
                                + (c1 * rand1[j] * (pbest_x[k].genes[j] - particles[k].genes[j]))
                                + (c2 * rand2[j] * (gbest_x.genes[j] - particles[k].genes[j]));
                        }

                        for j in 0..dim {
                            if v[k][j] > v_max[j] {
                                //index1.push(j);
                                v[k][j] = v_max[j];
                            }

                            if v[k][j] < v_min[j] {
                                //index2.push(j);
                                v[k][j] = v_min[j];
                            }
                        }

                        // Update particles positions
                        for j in 0..dim {
                            particles[k].genes[j] = particles[k].genes[j] + v[k][j];
                        }

                        for j in 0..dim {
                            if particles[k].genes[j] > ub[j] {
                                particles[k].genes[j] = ub[j];
                            }

                            if particles[k].genes[j] < lb[j] {
                                particles[k].genes[j] = lb[j];
                            }
                        }
                    }
                    cgcurve[t] = gbest_x.fitness.unwrap();
                }

                //return results
                let duration = chronos.elapsed();
                let result = OptimizationResult {
                    best_fitness: gbest_x.fitness,
                    best_genome: Some(gbest_x),
                    convergence_trend: Some(cgcurve),
                    computation_time: Some(duration),
                    err_report: None,
                };
                // copy result to PSO struct
                result
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct PSOparams<'a> {
    pub population_size: usize,
    pub problem_dimension: usize,
    pub max_iterations: usize,
    pub lower_bounds: &'a [f64],
    pub upper_bounds: &'a [f64],
    pub c1: f64,
    pub c2: f64,
}

#[allow(dead_code)]
impl<'a> PSOparams<'a> {
    pub fn new(
        p_size: usize,
        dim: usize,
        max_iter: usize,
        lb: &'a [f64],
        ub: &'a [f64],
        c1: f64,
        c2: f64,
    ) -> Result<PSOparams<'a>, String> {
        let params = PSOparams {
            population_size: p_size,
            problem_dimension: dim,
            max_iterations: max_iter,
            lower_bounds: lb,
            upper_bounds: ub,
            c1,
            c2,
        };

        match params.check() {
            Err(error) => Err(error),
            Ok(()) => Ok(params),
        }
    }
}

impl<'a> Parameters for PSOparams<'a> {
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

impl<'a> Default for PSOparams<'a> {
    ///
    /// Return default values of parameters, as following :
    ///
    /// ~~~
    ///
    /// use sefar::algos::pso::PSOparams;
    ///
    ///  PSOparams{
    ///     population_size : 10,
    ///     dimensions : 3,
    ///     max_iterations : 100,
    ///     lower_bounds : &[100.0f64, 100.0, 100.0],
    ///     upper_bounds : &[-100.0f64, -100.0, -100.0],
    ///     c1 : 2.0f64,
    ///     c2 : 1.0f64,
    /// };
    /// ~~~
    ///
    fn default() -> Self {
        PSOparams {
            population_size: 10,
            problem_dimension: 3,
            max_iterations: 100,
            lower_bounds: &[-100.0f64, -100.0, -100.0],
            upper_bounds: &[100.0f64, 100.0, 100.0],
            c1: 2.0f64,
            c2: 1.0f64,
        }
    }
}

impl<'a> Display for PSOparams<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Pop.Size: {}, Problem dim.: {}, Max.Iter: {}, c1: {}, c2: {}",
            self.population_size, self.problem_dimension, self.max_iterations, self.c1, self.c2
        )
    }
}

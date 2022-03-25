//
// Implementation of Particle Swarm Optimization algorithm (PSO)
// 
 
extern crate rand;
use rand::distributions::{Distribution, Uniform};
//use rand::prelude::ThreadRng;
use std::time::Instant;

use crate::core::eoa::EOA;
use crate::core::genome::Genome;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
use crate::core::optimization_result::OptimizationResult;
use crate::common::*;

pub struct PSO<'a, T: Problem> {
    problem : &'a mut T,
     
}




#[derive(Debug, Clone)]
pub struct PSOparams<'a> {
    pub population_size : usize,
    pub dimensions : usize,
    pub max_iterations: usize,
    pub lower_bounds : &'a [f64],
    pub upper_bounds : &'a [f64],
    pub c1 : f64,
    pub c2 : f64,
}

#[allow(dead_code)]
impl<'a> PSOparams<'a>{
    pub fn new(p_size: usize, dim : usize, max_iter : usize, lb : &'a [f64], 
    ub : &'a [f64], c1 :f64, c2 :f64)->PSOparams<'a> {
        PSOparams{
            population_size : p_size,
            dimensions : dim,
            max_iterations : max_iter,
            lower_bounds : lb,
            upper_bounds : ub,
            c1,
            c2,            
        }
    }
    
    ///
    /// Return default values of parameters, as following :
    /// 
    /// ~~~
    /// 
    ///  use sefar::sequential_algos::pso::*;
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
    pub fn default()->Self{
        PSOparams{
            population_size : 10,
            dimensions : 3,
            max_iterations : 100,
            lower_bounds : &[-100.0f64, -100.0, -100.0],
            upper_bounds : &[100.0f64, 100.0, 100.0],
            c1 : 2.0f64,
            c2 : 1.0f64,
        }
    }
}

impl<'a> Parameters for PSOparams<'a> {

   fn get_population_size(&self)->usize{
        self.population_size
    }

   fn get_dimensions(&self)-> usize {
        self.dimensions
    }

    fn get_max_iterations(&self)->usize{
        self.max_iterations
    }

    fn get_lower_bounds(&self)-> Vec<f64>{
        self.lower_bounds.to_vec()
    }

    fn get_upper_bounds(&self)-> Vec<f64>{
        self.upper_bounds.to_vec()
    }        
}

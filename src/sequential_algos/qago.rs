
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

///
/// QAGO
/// Reference:
/// 
/// 
/// 
/// 
/// 

#[derive(Debug)]
pub struct QAGO<'a, T : Problem> {
     pub problem : &'a mut T,
     pub params : &'a QAGOparams<'a>,
     pub optimization_result : OptimizationResult,
}

impl <'a, T : Problem> EOA for QAGO<'a, T>{
    fn run(&mut self)-> OptimizationResult {




        




        let result = OptimizationResult{
            best_genome : None,
            best_fitness : None, 
            convergence_trend : None,
            computation_time : None,
            err_report : None,
        };
        return result;  


    }
}






#[derive(Debug, Clone)]
pub struct QAGOparams<'a> {
    pub population_size : usize,
    pub dimensions : usize,
    pub max_iterations: usize,
    pub lower_bounds : &'a [f64],
    pub upper_bounds : &'a [f64],
}
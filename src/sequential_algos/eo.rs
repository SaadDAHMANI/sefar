//
// Implementation of Equilibrium Optimizer
// 
 
extern crate rand;
use rand::distributions::{Distribution, Uniform};
use std::time::Instant;

use crate::core::eoa::EOA;
use crate::core::genome::Genome;
use crate::sequential_algos::eo_params::EOparams;
use crate::core::problem::Problem;
use crate::core::optimization_result::OptimizationResult;
use crate::common::*;

///
/// Sequential Equilibrium Optimizer (EO)
/// 
#[derive(Debug)]
pub struct EO<'a, T : Problem> {
     pub problem : &'a mut T,
     pub params : &'a EOparams<'a>,
     pub optimization_result : OptimizationResult,
}

impl<'a, T : Problem> EO<'a, T>{

    pub fn new(settings :&'a EOparams, problem : &'a mut T )->Self{
       
        let result = OptimizationResult{
            best_genome : None,
            best_fitness :None,
            convergence_trend : None,
            computation_time : None,
            err_report : None, 
        };

        EO { 
             problem,
             params: settings,
             optimization_result: result, 
           }
    }
}

impl<'a, T: Problem> EOA for EO<'a, T> {
    
    fn run(&mut self)-> OptimizationResult{






        let result = OptimizationResult{
            best_genome : None,
            best_fitness :None,
            convergence_trend : None,
            computation_time : None,
            err_report : None, 
       };
       return result; 
        
    }

    fn run_epoch(&mut self){

    }

}


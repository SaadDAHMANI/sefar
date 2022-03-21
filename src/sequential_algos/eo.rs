//
// Implementation of Equilibrium Optimizer
// 
 
extern crate rand;
use rand::distributions::{Distribution, Uniform};
use rand::prelude::ThreadRng;
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

     // private 
     ceq1 : Vec<f64>,
     ceq2 : Vec<f64>,
     ceq3 : Vec<f64>,
     ceq4 : Vec<f64>,
     ceq_ave : Vec<f64>,
     ceq1_fit : f64,
     ceq2_fit : f64,
     ceq3_fit : f64,
     ceq4_fit : f64,

     // to store agents fitness values

     fitness : Vec<f64>,
     fit_old : Vec<f64>,
     c_old : Vec<Vec<f64>>,
     c_pool : Vec<Vec<f64>>,
     lambda : Vec<f64>,
     ceq : Vec<f64>, 
     r : Vec<f64>,
     r1 : Vec<f64>,
     r2 : Vec<f64>,
     f : Vec<f64>,
     
     interval : Uniform<usize>,
     rng : ThreadRng,
     c : Vec<Genome>,   

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

        let dim = settings.dimensions.clone();
        let particles_no = settings.population_size.clone();

        let _ceq1 = vec![0.0f64; dim];
        let _ceq2 = vec![0.0f64; dim];
        let _ceq3 = vec![0.0f64; dim];
        let _ceq4 = vec![0.0f64; dim];

        let _ceq_ave = vec![0.0f64; dim];
            
        let _ceq1_fit = f64::MAX;
        let _ceq2_fit = f64::MAX;
        let _ceq3_fit = f64::MAX;
        let _ceq4_fit = f64::MAX;

        let _fitness = vec![0.0f64; particles_no];
        let _fit_old = vec![0.0f64; particles_no];
         
        let _c_old = vec![vec![0.0f64; dim]; particles_no];
        let _c_pool = vec![vec![0.0f64; dim]; 5];
        let _lambda = vec![0.0f64; dim];
        let _ceq = vec![0.0f64; dim];

        let _r = vec![0.0f64; dim];
        let _r1 = vec![0.0f64; dim];
        let _r2 = vec![0.0f64; dim];
        let  _f = vec![0.0f64; dim];

        let _interval : Uniform<usize> = Uniform::from(0.._c_pool.len());
        let mut _rng = rand::thread_rng();

        let _c = initialize(&settings.clone());

        EO { 
             problem,
             params: settings,
             optimization_result: result, 
            
             // private 
             ceq1 : _ceq1,
             ceq2 : _ceq2,
             ceq3 : _ceq3,
             ceq4 : _ceq4,
             ceq_ave : _ceq_ave,
             ceq1_fit : _ceq1_fit,
             ceq2_fit : _ceq2_fit,
             ceq3_fit : _ceq3_fit,
             ceq4_fit : _ceq4_fit,
             fitness :_fitness,
             fit_old : _fit_old,
             c_old : _c_old,
             c_pool : _c_pool,
             lambda :_lambda,
             ceq :_ceq,
             r : _r,
             r1 : _r1,
             r2 : _r2,
             f :_f,
             interval : _interval,
             rng : _rng,
             c : _c,   
            }
    }

   

}

impl<'a, T: Problem> EOA for EO<'a, T> {
    
    fn run(&mut self)-> OptimizationResult{
        
        let mut iter : usize = 0;

        while iter < self.params.max_iterations {


             iter +=1;  
        }

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


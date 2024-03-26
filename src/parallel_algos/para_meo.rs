pub fn para_meo(){


}

//
// Implementation of Modified Equilibrium Optimizer (m-EO)
// 
 
extern crate rand;
//use rand::distributions::{Distribution, Uniform};

//use rand::prelude::ThreadRng;
//use std::time::Instant;

//use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};

use crate::core::eoa::{EOA, InitializationMode};
//use crate::core::genome::Genome;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
use crate::core::optimization_result::OptimizationResult;
//use crate::common::*;
use crate::sequential_algos::eo::EOparams;

///
/// Sequential Modified Equilibrium Optimizer (MEO)
/// Reference:
/// "Gupta, S., Deep, K., & Mirjalili, S. (2020).
/// An efficient equilibrium optimizer with mutation strategy for numerical optimization. 
/// Applied Soft Computing, 96, 106542."
/// 
#[derive(Debug)]
pub struct ParaMEO<'a, T : Problem> {
     pub problem : &'a mut T,
     pub params : &'a EOparams<'a>,
     pub optimization_result : OptimizationResult,
}

impl<'a, T : Problem> ParaMEO<'a, T>{

    pub fn new(settings :&'a EOparams, problem : &'a mut T )->Self{       
        let result = OptimizationResult{
            best_genome : None,
            best_fitness :None,
            convergence_trend : None,
            computation_time : None,
            err_report : None, 
        };
       
        ParaMEO { 
             problem,
             params: settings,
             optimization_result: result,            
        }
    }   
}


impl<'a, T: Problem> EOA for ParaMEO<'a, T> {
   
    fn run(&mut self)-> OptimizationResult{
        
        //let chronos = Instant::now();
    
        //check paramaters
        //let params = self.params.clone();

        //__________________________________________________________________

    

      
        //___________________________________________________________________    

        
        match self.params.check(){
            Err(error) => OptimizationResult::get_none(error),
            Ok(()) =>{

                 
                // Step 1: initialize the population randomly within the solution space 
                let c = self.initialize(self.params, InitializationMode::RealUniform);
                
                // Step 2 : Evaluate the fitness value of each candidate soluion
                //for genom in c.iter_mut(){
                
                // c.par_iter_mut().enumerate().for_each(|(_i,g)| g.fitness =
                // Some(self.problem.objectivefunction(&g.genes)));
        
                
                for genom in c.iter(){
                    println!("id : {}, fitness : {:?}", genom.id, genom.fitness);
                }
           


                //return results                
                let result = OptimizationResult{
                    best_genome : None,
                    best_fitness : None,
                    convergence_trend : None,
                    computation_time : None,
                    err_report : None,
                };
                // copy result to EO struct
                self.optimization_result = result.clone();
                result
            }
        }
    }    
}    
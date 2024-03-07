/// -------------------------------------------------------------------------------------
/// Plant competition optimization (PCO)
/// Implemented in Rust programming language by Saad Dahmani (s.dahmani@univ-bouira.dz)
/// Original Matlab code :
/// https://github.com/iman-aliabdi/PCO-Plant-Competition-Optimization
/// Plant competition optimization: A novel metaheuristic algorithm
/// Piper's link : https://onlinelibrary.wiley.com/doi/abs/10.1111/exsy.12956
/// -------------------------------------------------------------------------------------
/// 

extern crate rand;
use rand::distributions::{Distribution, Uniform};

use crate::core::eoa::EOA;
use crate::core::genome::Genome;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
use crate::core::optimization_result::OptimizationResult;
use crate::common::*;

///
/// Plant competition optimization (PCO) (Sequential)
/// Reference:
/// Plant competition optimization: A novel metaheuristic algorithm
/// Piper's link : https://onlinelibrary.wiley.com/doi/abs/10.1111/exsy.12956
/// Original Matlab code:
/// https://github.com/iman-aliabdi/PCO-Plant-Competition-Optimization
/// 
#[derive(Debug)]
pub struct PCO<'a, T : Problem> {
     pub problem : &'a mut T,
     pub params : &'a PCOparams<'a>,
     pub optimization_result : OptimizationResult,
}

impl<'a, T : Problem> PCO<'a, T>{

    pub fn new(settings :&'a PCOparams, problem : &'a mut T )->Self{       
        let result = OptimizationResult{
            best_genome : None,
            best_fitness :None,
            convergence_trend : None,
            computation_time : None,
            err_report : None, 
        };
       
        PCO { 
             problem,
             params: settings,
             optimization_result: result,            
        }
    }   
}


impl<'a, T: Problem> EOA for PCO<'a, T> {
   
    fn run(&mut self)-> OptimizationResult{


        let dim : usize = self.params.dimensions;
        let ub = self.params.upper_bounds;
        let lb = self.params.lower_bounds;

        //----------------------------------------
        let mut rmax : Vec<f64> = vec![0.0f64; dim];




        //-----------------------------------------
        
        let mut i : usize = 0;
        for (a,b) in ub.iter().zip(lb.iter()){
            rmax[i] = a-b;
            i+=1;
        }








        let result = OptimizationResult{
            best_genome : None, //Some(Genome::new(0, self.params.dimensions)),
            best_fitness : None, //Some(-1111.1111),
            convergence_trend : None, //Some(convergence_curve),
            computation_time : None, //Some(duration),
            err_report : None,
        };
        return result;   

    }    
}







#[derive(Debug, Clone)]
pub struct PCOparams<'a> {
    /// Number of initial plants.
    pub population_size : usize,

    /// Size of the problem (number of decision variables).
    pub dimensions : usize,

    /// Maximum number of iterations.
    pub max_iterations: usize,

    pub lower_bounds : &'a [f64],
    pub upper_bounds : &'a [f64],
    
    /// vmax : Maximum size of plants.
    pub vmax : usize,

    /// Maximum of Plants Number.
    pub max_plant_number : usize,

    /// Growth rate.
    pub theta : f64,
    
    /// Growth parameter. 
    pub k : f64,

    ///Seed migration rate.
    pub miu : f64,    
}

#[allow(dead_code)]
impl<'a> PCOparams<'a>{
    pub fn new(p_size: usize, dim : usize, max_iter : usize, lb : &'a [f64], 
    ub : &'a [f64], vmax : usize, max_plant_number : usize, theta : f64, k : f64, miu : f64)-> Result<PCOparams<'a>, String> {
                      
        let params = PCOparams{
            population_size : p_size,
            dimensions : dim,
            max_iterations : max_iter,
            lower_bounds : lb,
            upper_bounds : ub,
            vmax,
            max_plant_number,
            theta,
            k,
            miu,
        };

       match params.check() {
           Err(error)=> Err(error),
           Ok(())=> Ok(params),
       }
    }
}

impl<'a> Parameters for PCOparams<'a> {

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
 

impl<'a> Default for PCOparams<'a>{

    ///
    /// Return default values of parameters, as following :
    /// 
    /// ~~~
    /// 
    ///  use sefar::sequential_algos::pco::*;
    /// 
    ///  PCOparams{
    ///     population_size : 10,
    ///     dimensions : 3,
    ///     max_iterations : 100,
    ///     lower_bounds : &[-100.0f64, -100.0, -100.0],
    ///     upper_bounds : &[100.0f64, 100.0, 100.0],
    ///     vmax : 20,
    ///     max_plant_number : 500,
    ///     theta : 0.005,
    ///     k : 0.1,
    ///     miu : 0.05,
    /// };
    /// ~~~
    /// 
    fn default()->Self{
        PCOparams{
            population_size : 10,
            dimensions : 3,
            max_iterations : 100,
            lower_bounds : &[-100.0f64, -100.0, -100.0],
            upper_bounds : &[100.0f64, 100.0, 100.0],
            vmax : 20,
            max_plant_number : 500,
            theta : 0.005,
            k : 0.1,
            miu : 0.05,
        }
    }
}


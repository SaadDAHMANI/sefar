use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::genome::Genome;
use rand::distributions::{Distribution, Uniform};

///
/// Public trait for Evolutionary Optimization Algorithms
/// 
pub trait EOA {
    
    fn initialize<P: Parameters>(&self, params: &P)-> Vec<Genome>{

       let n: usize = params.get_population_size();
       let dim: usize = params.get_dimensions();
       let lb = params.get_lower_bounds();
       let ub = params.get_upper_bounds();
    
       let mut positions : Vec<Genome> = Vec::with_capacity(n);

       let intervall01 = Uniform::from(0.0f64..=1.0f64);
       let mut rng = rand::thread_rng();              
  
       for i in 0..n{
             let mut sln = Genome::new(i, dim); 
            
             for  j in 0..dim {   
              //  positions[i][j]= intervall01.sample(&mut rng)*(ub-lb)+lb;                         
              sln.genes[j]= intervall01.sample(&mut rng)*(ub[j]-lb[j])+lb[j];   
            }
            positions.push(sln);
        }        
        positions
    }    
    ///
    /// Run algorithm until reach stopping criterion
    /// 
    fn run(&mut self)-> OptimizationResult {
        OptimizationResult {
             best_genome : None,
             best_fitness :None,
             convergence_trend : None,
             computation_time : None,
             err_report : None, 
        } 
    }  
        
    fn check_parameters<P: Parameters>(params : &P)-> Result<(), String> {
    
        let mut errors : usize = 0;
        let mut msg : String = String::new();
        
        if params.get_population_size() == 0 {
            msg = String::from("population_size must be greater than 0!; \n");
            errors +=1;
        }
    
        if params.get_dimensions() ==0 {
            msg = format!("{} Search space dimensions must be greater than 0!; \n", msg);
            errors +=1;
        }
    
        if params.get_max_iterations() ==0 {
            msg = format!("{} Iterations count (max_iterations) must be greater than 0!; \n", msg);
            errors +=1;
        }
    
        if params.get_lower_bounds().is_empty() {
            msg = format!("{} Lower_bounds length must be greater than 0!; \n", msg);
            errors +=1;
        }
    
        if params.get_upper_bounds().is_empty() {
            msg = format!("{} Upper_bounds length must be greater than 0!; \n", msg);
            errors +=1;
        }
    
        if params.get_lower_bounds().len() != params.get_upper_bounds().len() {
            msg = format!("{} Lower_bounds & Upper_bounds lengths must be equal!; \n", msg);
            errors +=1;
        }
    
        if params.get_lower_bounds().len() != params.get_dimensions() || params.get_upper_bounds().len() != params.get_dimensions() {
            msg = format!("{} Lower_bounds & Upper_bounds lengths must equal search space dimension!; \n", msg);
            errors +=1;
        }
       
        if errors > 0  {
            msg = format!("There are [{}] errors : \n {}", errors, msg.trim());       
            Err(msg)
        }
        else {
            Ok(())
        }      
    }   


}
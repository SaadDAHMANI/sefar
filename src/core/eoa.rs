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
       let lb = &params.get_lower_bounds();
       let ub = &params.get_upper_bounds();
    
       let mut positions : Vec<Genome> = Vec::with_capacity(n);

       let intervall01 = Uniform::from(0.0f64..=1.0f64);
       let mut rng = rand::thread_rng();              
  
       for i in 0..n{
             let mut sln = Genome::new(i, dim); 
            
             for  j in 0..dim {   
              //  positions[i][j]= intervall01.sample(&mut rng)*(ub-lb)+lb;                         
              sln.genes[j]= intervall01.sample(&mut rng)*(ub[j]-lb[j]) + lb[j];   
            }
            positions.push(sln);
        }        
        positions
    }    
    
    ///
    /// Run algorithm until reach stopping criterion and return optiization result
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
        
   
    fn randomize(randvect : &mut Vec<f64>) {    
        let between = Uniform::from(0.0..=1.0);
        let mut rng = rand::thread_rng();
                
        for item in randvect.iter_mut() {
            *item = between.sample(&mut rng);
        }     
    }


}
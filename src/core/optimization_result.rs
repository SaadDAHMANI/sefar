//use std::fmt::Display;
use std::time::Duration;
use crate::core::genome::Genome;

#[derive(Debug, Clone)]
pub struct OptimizationResult {
    pub best_genome : Option<Genome>,
    pub best_fitness : Option<f64>,
    pub convergence_trend : Option<Vec<f64>>,
    pub computation_time : Option<Duration>,
    pub err_report : Option<String>,         
 }
 
 #[allow(dead_code)]
 impl OptimizationResult{
    
    pub fn get_none(msg : String)->OptimizationResult {
        OptimizationResult{
            best_genome : None,
            best_fitness : None,
            convergence_trend : None,
            computation_time : None,
            err_report : Some(msg),                              
        }   
     }
 }  

/*  impl Display for OptimizationResult{
     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
         write!(f, "Best-fitness : {:?}; Best-solution : {:?}; Time : {:?}; Err-report: {:?}", self.best_fitness, self.best_genome, self.computation_time, self.err_report)
     }
 } */

 impl ToString for OptimizationResult{
    fn to_string(&self) -> String {
        format!("Best-fitness : {:?} \n; Best-solution : {:?} \n; Time : {:?} \n; Err-report: {:?}", self.best_fitness, self.best_genome, self.computation_time, self.err_report)
    }
 }
 

 
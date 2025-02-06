use std::path::{self, Path};
//use std::fmt::Display;
use std::time::Duration;
use std::io::{Read, Write};
use std::fs::{exists, File};

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

 impl OptimizationResult{
    pub fn save(&self, header :Option<&str>, filename : &str)->std::io::Result<()> {
      
        let mut file : File = match std::fs::metadata(filename) {
            Err(_eror) =>{
                let file = File::create(&filename)?;
                file
            },

            Ok(mdata) =>{
                if mdata.is_file(){
                    let  file : File = File::open(&filename)?;
                    file
                }
                else{ let file = File::create(&filename)?;
                    file
                 }
            }
        };

        match header {
            None =>{},
            Some(header)=> {writeln!(file, "{:?}", header)?;},
        };
        
        writeln!(file, "Best_fitness : {:?}", self.best_fitness)?;        

        Ok(())

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
 

 
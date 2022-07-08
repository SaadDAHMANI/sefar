
use sefar::core::problem::Problem;
//use crate::core::genome::Genome;

#[derive(Debug, Clone)]
pub struct Regression {
   pub learn_in : Vec<Vec<f64>>,
   pub learn_out : Vec<f64>,
   pub test_in : Vec<Vec<f64>>,
   pub test_out : Vec<Vec<f64>>,
   pub file : String,
}
impl Regression {
    pub fn new(fileptah : String)-> Regression {
        
        //let root = String::from("/home/sd/Documents/Rust_apps/sefar/src/bin/data");
        
        //let path = format!("{}/{}", root, "Coxs_data.csv"); 
     
        let mut incols = Vec::new();
    
        incols.push(2); 
        incols.push(3);  
        //incols.push(4);  
        //incols.push(5);
        //incols.push(6); 
        //incols.push(7);
        //incols.push(8);
        //incols.push(9);
        //incols.push(10);
        //incols.push(11);
        //incols.push(12); 
        //incols.push(13);
        //incols.push(14); 
        //incols.push(15);
        //incols.push(16); 
        //incols.push(17);
        //incols.push(18);
        //incols.push(19);

    let mut outcols = Vec::new();
    outcols.push(1);
    
    let learn_part : usize = 4233; // Total =6047; Test = 1814 (70%, 30%) 
    
    let ds0 =  Dataset::read_from_csvfile(&fileptah, &incols, &outcols);

    let (ds_learn, ds_test) = ds0.split_on_2(learn_part);

    let learn_in = ds_learn.inputs;
    let learn_out = Dataset::get_first_items(&ds_learn.outputs);
    let test_in = ds_test.inputs;
    let test_out = ds_test.outputs;    
        
    Regression{
        learn_in,
        learn_out,
        test_in,
        test_out,
        file : fileptah,
        }
    }   
}

impl Problem for Regression {
    fn objectivefunction(&mut self, genome : &[f64]) ->f64 {
              
        let n = self.learn_in.len();
        let mut computed : Vec<f64> = Vec::with_capacity(n);
        
        let m = genome.len()-1;

        for i in 0..n {
            
            let mut q =0.0f64;

            for j in 0..m {
                 q+= genome[j]*self.learn_in[i][j];     
            }
            q += genome[m];

            computed[i] = q;
        }
     
        let fit = Dataset::compute_rmse(&computed, &self.learn_out);
        
        let fitness = match fit{
            None => f64::powi(10.0f64,10),
            Some(fit)=> fit,
        };

        fitness
    }
}

use sefar::core::problem::Problem;
//use crate::core::genome::Genome;

#[derive(Debug, Clone)]
pub struct Regression {
   pub learn_in : Vec<Vec<f64>>,
   pub learn_out : Vec<f64>,
   pub test_in : Vec<Vec<f64>>,
   pub test_out : Vec<f64>,
   pub file : String,
}
impl Regression {
    pub fn new(fileptah : String)-> Regression {      
            
    let incols = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]; // M10

    let incols = incols.to_vec();
     
    let mut outcols = Vec::new();
    outcols.push(1);
    
    let learn_part : usize = 4233; // Total =6047; Test = 1814 (70%, 30%) 
    
    let ds0 =  Dataset::read_from_csvfile(&fileptah, &incols, &outcols);

    let (ds_learn, ds_test) = ds0.split_on_2(learn_part);
     
    //println!("ds_learn.outputs = {:?}", ds_learn.outputs);

    let learn_in = ds_learn.inputs;
    let learn_out = Dataset::get_first_items(&ds_learn.outputs);

    //println!("learn_out = {:?}", learn_out);

    let test_in = ds_test.inputs;
    let test_out = Dataset::get_first_items(&ds_test.outputs);    
        
    Regression{
        learn_in,
        learn_out,
        test_in,
        test_out,
        file : fileptah,
        }
    }   

    pub fn compute_learn_indexes(&self, params : &Vec<f64>)->(f64, f64) {

        let n = self.learn_in.len();
        let mut computed : Vec<f64> = Vec::new();
        
        let m = self.learn_in[0].len();
        
        //println!("m= {:?}", m);

        for i in 0..n {            
            let mut q =0.0f64;

            for j in 0..m {
                 q+= params[j]*self.learn_in[i][j];     
            }
            q += params[m];

            computed.push(q);
        }
     
        let rmse = Dataset::compute_rmse(&computed, &self.learn_out);
        
        let rmsel = match rmse{
            None => -1.0f64,
            Some(value)=> value,
        };

        let r2 = Dataset::compute_determination_r2(&computed, &self.learn_out);

        let r2l = match r2 {
            None => -1.0f64,
            Some(value)=> value,
        };

        (rmsel, r2l)
    }

    pub fn compute_test_indexes(&self, params : &Vec<f64>)->(f64, f64) {

        let n = self.test_in.len();
        let mut computed : Vec<f64> = Vec::new();
        
        let m = self.test_in[0].len();
        
        //println!("m= {:?}", m);

        for i in 0..n {            
            let mut q =0.0f64;

            for j in 0..m {
                 q+= params[j]*self.test_in[i][j];     
            }
            q += params[m];

            computed.push(q);
        }
     
        let rmse = Dataset::compute_rmse(&computed, &self.test_out);
        
        let rmset = match rmse{
            None => -1.0f64,
            Some(value)=> value,
        };

        let r2 = Dataset::compute_determination_r2(&computed, &self.test_out);

        let r2t = match r2 {
            None => -1.0f64,
            Some(value)=> value,
        };

        (rmset, r2t)
    }
    
    pub fn compute_result_indexes(&self, params : &Vec<f64>)->(f64, f64, f64, f64) {
        let (rmsel, r2l)=self.compute_learn_indexes(params);
        let (rmset, r2t)=self.compute_test_indexes(params);
        (rmsel, rmset, r2l, r2t)
    }


}

impl Problem for Regression {
    fn objectivefunction(&mut self, genome : &[f64]) ->f64 {
              
        let n = self.learn_in.len();
        let mut computed : Vec<f64> = Vec::new();
        
        let m = self.learn_in[0].len();
        
        //println!("m= {:?}", m);

        for i in 0..n {            
            let mut q =0.0f64;

            for j in 0..m {
                 q+= genome[j]*self.learn_in[i][j];     
            }
            q += genome[m];

            computed.push(q);
        }
     
        let fit = Dataset::compute_rmse(&computed, &self.learn_out);
        
        let fitness = match fit{
            None => f64::powi(10.0f64,10),
            Some(fit)=> fit,
        };

        fitness
    }
}


use std::error::Error;

use csv;

use csv::Writer;

//use serde::de::DeserializeOwned;
//use serde::Deserialize;

#[derive(Debug, Clone)]
pub struct Dataset {
    pub inputs : Vec<Vec<f64>>,
    pub outputs : Vec<Vec<f64>>, 
    pub file_path : Option<String>,    
    //inputs_headers : &'a mut Vec<String>,
    //outputs_headers : &'a mut Vec<String>,
}

impl Dataset{

    pub fn new(inputs : Vec<Vec<f64>>, outputs :Vec<Vec<f64>>, file : Option<String>)-> Dataset {
        Dataset{
            inputs : inputs,
            outputs : outputs, 
            file_path : file,            
        }
    }
    
    pub fn len(&self)->(usize, usize){
      	(self.inputs.len(), self.outputs.len())
    }
    
    pub fn clone_data(&self)->(Vec<Vec<f64>>, Vec<Vec<f64>>){
    	  (self.inputs.clone(), self.outputs.clone())    	
    }
    
    pub fn read_from_csvfile(path : &String, in_cols : &Vec<usize>, out_cols : &Vec<usize>)->Dataset{
            
        let origing_data = Dataset::readall_from_file(path);

        let data = match origing_data {
            Ok(data) => data,
            Err(error) => panic!("Problem was found : {:?}", error),
        };

        let mut input = Vec::new();
        let mut output = Vec::new();

        for i in 0..data.len() {
            let mut row_in = Vec::new();
            let mut row_out = Vec::new();

            for j in in_cols.iter(){
                if *j< data[i].len() {
                    row_in.push(data[i][*j]);
                }  
            }

            for k in out_cols.iter() {
                if *k < data[i].len() {
                    row_out.push(data[i][*k]);
                }
            }
 
            input.push(row_in);
            output.push(row_out);
        }
        
        let flepath = path.clone();
        let ds = Dataset {
            inputs : input,
            outputs : output,
            file_path : Some(flepath),
        };
        return ds;
    }

    pub fn readall_from_file(path : &str)-> Result< Vec<Vec<f64>>, Box<dyn Error>> {
    
         let mut reader = csv::Reader::from_path(path)?;
    
         let headers = reader.headers()?;
    
         println!("Headers :  {:?}", headers);

         let mut data = Vec::new();

        for result in reader.records() {

             let record = result?;

             if record.is_empty()==false {
                 let l = record.len();
                 let mut row = vec![0.0f64; l];

                 for i in 0..l {
                 row[i] = record[i].parse()?; 
                }
            data.push(row)            
            }       
        }
         Ok(data)
    }
   
    pub fn write_to_csv(path : &String, header : &Option<String>, data : &Vec<f64>)-> Result<(), Box<dyn Error> > {
        let mut wtr = Writer::from_path(path)?;

        match header {
            Some(header) => wtr.write_record(&[header])?,
            None => (),
        };

        for i in 0..data.len() {
             wtr.write_record(&[data[i].to_string()])?;
            
        }

        wtr.flush()?;
        
        Ok(())   
    }

    pub fn write_to_csv2(path : &String, headers : &Option<Vec<String>>, data : &Vec<Vec<f64>>)-> Result<(), Box<dyn Error> > {
        let mut wtr = Writer::from_path(path)?;

        match headers {
            Some(header) => wtr.write_record(header.iter())?,
            None => (),
        };

        for row in data {
             let cols_str: Vec<_> = row.iter().map(ToString::to_string).collect();   
             //let line = cols_str.join(",");
             wtr.write_record(cols_str.iter())?;
        }

        wtr.flush()?;
        
        Ok(())   
    }

    ///
    /// shuffled data in the intervalle [0, 1] 
    /// 
    pub fn get_shuffled(&self)->Dataset {
        let mut maxin = Vec::new();

        if self.inputs.len()> 0 {
             maxin = self.inputs[0].clone();        
        }        

        let mut maxout = Vec::new();
        if self.outputs.len()> 0 {
             maxout = self.outputs[0].clone();
        }
        
        let icountin = self.inputs.len();
        let jcountin = self.inputs[0].len(); 

        for j in 0.. jcountin {
            for i in 0..icountin {
               if maxin[j] < self.inputs[i][j] {
                   maxin[j]=self.inputs[i][j];
               }  
            }   
        }

        let icountout = self.outputs.len();
        let jcountout = self.outputs[0].len(); 

        for j in 0.. jcountout {
            for i in 0..icountout {
               if maxout[j] < self.outputs[i][j] {
                   maxout[j]=self.outputs[i][j];
               }  
            }   
        }

        let mut shuffled = self.clone();
        for i in 0..icountin {
            for j in 0.. jcountin {
                if maxin[j] != 0.0 {
                    shuffled.inputs[i][j]= shuffled.inputs[i][j]/maxin[j]; 
                }                               
            }
        }
        
        for i in 0..icountout {
            for j in 0.. jcountout {
                if maxout[j] != 0.0 {
                    shuffled.outputs[i][j]= shuffled.outputs[i][j]/maxout[j]; 
                }                               
            }
        }

        return shuffled;
    } 

     ///
    /// shuffled data in the intervalle [0, 0.9] 
    /// 
    pub fn get_shuffled_09(&self)->Dataset {
        let mut maxin = Vec::new();
        let mut minin = Vec::new();

        if self.inputs.len()> 0 {
             maxin = self.inputs[0].clone();
             minin =self.inputs[0].clone();        
        }        

        let mut maxout = Vec::new();
        let mut minout = Vec::new();
        if self.outputs.len()> 0 {
             maxout = self.outputs[0].clone();
             minout = self.outputs[0].clone();
        }
        
        let icountin = self.inputs.len();
        let jcountin = self.inputs[0].len(); 

        // search min and max values
        for j in 0.. jcountin {
            for i in 0..icountin {
               if maxin[j] < self.inputs[i][j] {
                   maxin[j] = self.inputs[i][j];
               }
               
               if minin[j] > self.inputs[i][j] {
                   minin[j] = self.inputs[i][j]
               }   
            }   
        }

        let icountout = self.outputs.len();
        let jcountout = self.outputs[0].len(); 

         // search min and max values
        for j in 0.. jcountout {
            for i in 0..icountout {
               if maxout[j] < self.outputs[i][j] {
                     maxout[j]=self.outputs[i][j];
               }  

               if minout[j] > self.outputs[i][j] {
                     minout[j]=self.outputs[i][j];
                }  
            }   
        }

        let mut shuffled = self.clone();
        for i in 0..icountin {
            for j in 0.. jcountin {
                if maxin[j] != 0.0 {
                    shuffled.inputs[i][j]= 0.9 *(shuffled.inputs[i][j]-minin[j])/(maxin[j]-minin[j]); 
                }                               
            }
        }
        
        for i in 0..icountout {
            for j in 0.. jcountout {
                if maxout[j] != 0.0 {
                    shuffled.outputs[i][j]= 0.9*(shuffled.outputs[i][j]-minout[j])/(maxout[j]-minout[j]); 
                }                               
            }
        }

        return shuffled;
    } 

    pub fn split_on_2(&self, firstelemntscount : usize)->(Dataset, Dataset) {
        let totalcount = self.inputs.len();

        let firstcount = usize::min(firstelemntscount, totalcount);
        
        let mut datain1 = Vec::new();
        let mut datain2 = Vec::new();
        let mut dataout1 = Vec::new();
        let mut dataout2 = Vec::new();

        for i in 0..firstcount {
            datain1.push(self.inputs[i].clone());
            dataout1.push(self.outputs[i].clone());
        }

        for i in firstcount..totalcount {
            datain2.push(self.inputs[i].clone());
            dataout2.push(self.outputs[i].clone());
        } 

        
        let ds1 = Dataset::new(datain1, dataout1, None);
        let ds2 = Dataset::new(datain2, dataout2, None);

        (ds1, ds2)
    }  

    pub fn compute_rmse(ds1 : &Vec<f64>, ds2 : &Vec<f64>)-> Option<f64> {
         let mut rmse : f64 =0.0;

         let count = usize::min(ds1.len(), ds2.len());

         if count > 0 {
             for i in 0..count {
                 rmse += f64::powi(ds1[i]- ds2[i], 2);
             }
             rmse = rmse /count as f64;
             rmse = f64::sqrt(rmse); 
             Some(rmse)
         }
         else {
             None
         }
         
    }

    pub fn compute_correlation_r(ds1 : &Vec<f64>, ds2 : &Vec<f64>)->Option<f64> {
        let mut sum1 : f64 = 0.0;
        let mut sum2 : f64 = 0.0;
        let mut av1 : f64 = 0.0;
        let mut av2 : f64 = 0.0;
        let mut numerator : f64 = 0.0;

        let count = ds1.len();
       
        if ds1.len() == ds2.len() {
            if count > 0 {
                for i in 0..count {
                    av1 += ds1[i];
                    av2 += ds2[i];            
                };
    
                av1 = av1/count as f64;
                av2 = av2/count as f64;
    
                for i in 0..count {
                    sum1 += f64::powi(ds1[i] - av1, 2);
                    sum2 += f64::powi(ds2[i] - av2, 2);
                    numerator += (ds1[i] - av1)*(ds2[i] - av2) ;                
                };
                let denomenator = f64::sqrt(sum1*sum2);
                 
                let rvalue = numerator/denomenator;
                
                Some(rvalue)
            }
            else {None}
        }
        else {
            None
        }
       
    }

    pub fn compute_determination_r2(ds1 : &Vec<f64>, ds2 : &Vec<f64>)->Option<f64>{
       let r2 = match Dataset::compute_correlation_r(&ds1, &ds2) {
            Some(r)=> Some(f64::powi(r,2)),
            None => None,
        };
        r2
    } 

}


#[cfg(test)]
mod tests {
    use super::*;
   
    #[test]
    fn compute_correlation_r_test1() {
          let mut ds1 = Vec::new();
          ds1.push(1.2f64);
          ds1.push(2.2f64);
          ds1.push(3.2f64);  
          ds1.push(2.2f64);
          ds1.push(10.2f64);  

          let ds2 = ds1.clone();
           
          assert_eq!(Dataset::compute_correlation_r(&ds1, &ds2), Some(1.0));
    }

    #[test]
    fn compute_correlation_r_test2() {
          let mut ds1 = Vec::new();
          ds1.push(1.2f64);
          ds1.push(2.2f64);
          ds1.push(3.2f64);  
          ds1.push(2.2f64);
          ds1.push(10.2f64);  

          let ds2 = ds1.clone();

          ds1[3]=4.2;
           
          assert_ne!(Dataset::compute_correlation_r(&ds1, &ds2), Some(1.0));
    }

    #[test]
    fn compute_correlation_r_test3() {
        let mut ds1 = Vec::new();
        ds1.push(1.0f64);
        ds1.push(2.0f64);
        ds1.push(3.0f64);  
        ds1.push(4.0f64);
        ds1.push(5.0f64); 
        ds1.push(6.0f64);   

        let ds2 = ds1.clone();

        ds1[0]=2.0;
         
        assert_eq!(Dataset::compute_correlation_r(&ds1, &ds2), Some(0.9819805060619656));
  }

  #[test]
  fn compute_correlation_r_test4() {
      let mut ds1 = Vec::new();
      ds1.push(1.0f64);
      ds1.push(2.0f64);
      ds1.push(3.0f64);  
      ds1.push(10.0f64);
      ds1.push(11.0f64); 
      ds1.push(12.0f64);   

      let mut ds2 = Vec::new();
      ds2.push(2.0f64);
      ds2.push(2.0f64);
      ds2.push(3.0f64);  
      ds2.push(4.0f64);
      ds2.push(5.0f64); 
      ds2.push(6.0f64);   
      
      assert_eq!(Dataset::compute_correlation_r(&ds1, &ds2), Some(0.9533961104526778));
}

    #[test]
    fn compute_correlation_rmse_test1() {
        let mut ds1 = Vec::new();
          ds1.push(1.2f64);
          ds1.push(2.2f64);
          ds1.push(3.2f64);  
          ds1.push(2.2f64);
          ds1.push(10.2f64);  

          let ds2 = ds1.clone();
           
          assert_eq!(Dataset::compute_rmse(&ds1, &ds2), Some(0.0));
    }

    #[test]
    fn compute_correlation_rmse_test2() {
        let mut ds1 = Vec::new();
          ds1.push(1.2f64);
          ds1.push(2.2f64);
          ds1.push(3.2f64);  
          ds1.push(2.2f64);
          ds1.push(10.2f64);  

          let ds2 = ds1.clone();

          ds1[3]=5.2;
           
          assert_ne!(Dataset::compute_rmse(&ds1, &ds2), Some(0.0));
    }

    #[test]
    fn compute_correlation_rmse_test3() {
        let mut ds1 = Vec::new();
        ds1.push(1.0f64);
        ds1.push(2.0f64);
        ds1.push(3.0f64);  
        ds1.push(4.0f64);
        ds1.push(5.0f64); 
        ds1.push(6.0f64);   

        let ds2 = ds1.clone();

        ds1[0]=2.0;
         
        assert_eq!(Dataset::compute_rmse(&ds1, &ds2), Some(0.408248290463863));
  }

  #[test]
  fn compute_rmse_test4() {
      let mut ds1 = Vec::new();
      ds1.push(1.0f64);
      ds1.push(2.0f64);
      ds1.push(3.0f64);  
      ds1.push(10.0f64);
      ds1.push(11.0f64); 
      ds1.push(12.0f64);   

      let mut ds2 = Vec::new();
      ds2.push(2.0f64);
      ds2.push(2.0f64);
      ds2.push(3.0f64);  
      ds2.push(4.0f64);
      ds2.push(5.0f64); 
      ds2.push(6.0f64);   
      
      assert_eq!(Dataset::compute_rmse(&ds1, &ds2), Some(4.262237284181474));
}   

}

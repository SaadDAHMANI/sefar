include!("strm_regression.rs");
include!("dataset.rs");


use sefar::core::eoa::EOA;
use sefar::benchmarks::functions::{Sphere, F2};
use sefar::core::optimization_result::OptimizationResult;
use sefar::sequential_algos::eo::{EO, EOparams};
use sefar::sequential_algos::pso::{PSO, PSOparams};
//use sefar::sequential_algos::pco::{PCO, PCOparams};

use sefar::sequential_algos::meo::MEO;
use sefar::sequential_algos::qago::{QAGOparams, QAGO};
use sefar::sequential_algos::go::{GOparams, GO};

const DIM : usize = 7;
const POP_SIZE : usize = 20;
const KMAX : usize = 500; //1000*DIM/POP_SIZE;

fn main() {
    println!("Hello, sefar !");

    println!("Evaluation with Max_Iter = {}", KMAX);
    println!("______________________GO : F1______________________");

    go_f1_test1();

    println!("______________________GO : F2______________________");
    go_f2_test1();

    println!("_______________________QAGO : F1______________________");

    //qago_f1_test1();

    //pco_f1_test1();
    
    //eo_f1_test1();

  

    //peo_f1_test1();

    //meo_test1();

    //do_regression();
                                             

  }


  #[allow(dead_code)]
  fn go_f1_test1(){
      let mut settings : GOparams = GOparams::default();
      
      settings.population_size = POP_SIZE;
      settings.dimensions = DIM ;    
      settings.max_iterations = KMAX; 
  
         
      let lb =vec![-100.0f64; DIM];
      let ub =vec![100.0f64; DIM];
  
      settings.lower_bounds = lb.as_slice();
      settings.upper_bounds = ub.as_slice();  
  
      let mut fo = Sphere{};
  
      let mut algo : GO<Sphere> = GO::new(&settings, &mut fo);
      
      let result : OptimizationResult = algo.run();
  
      /*match result.convergence_trend{
          None => println!("QAGO: no convergence trend !!!"),
          Some(cv) => println!("QAGO: Convergence trend :\n {:?}", cv),
      };
      */
  
      /*match result.best_genome {
          None => println!("QAGO: no best solution !"),
          Some(bg)=> println!("QAGO: best-genome {:?}", bg),
      };
      */
  
      println!("Growth optimizer (GO) : F1 (Sphere) test; Result: {:?}", result.to_string());
}


#[allow(dead_code)]
fn go_f2_test1(){
    let mut settings : GOparams = GOparams::default();
    
    settings.population_size = POP_SIZE;
    settings.dimensions = DIM ;    
    settings.max_iterations = KMAX; 
       
    let lb =vec![-100.0f64; DIM];
    let ub =vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();  

    let mut fo = F2{};

    let mut algo : GO<F2> = GO::new(&settings, &mut fo);
    
    let result = algo.run();

    /*match result.convergence_trend{
        None => println!("QAGO: no convergence trend !!!"),
        Some(cv) => println!("QAGO: Convergence trend :\n {:?}", cv),
    };
    */

    /*match result.best_genome {
        None => println!("QAGO: no best solution !"),
        Some(bg)=> println!("QAGO: best-genome {:?}", bg),
    };
    */

    println!("Growth Optimizer (GO) : F2 test; Result: {:?}", result.to_string());
}
  


#[allow(dead_code)]
fn qago_f1_test1(){
    let mut settings : QAGOparams = QAGOparams::default();
    
    settings.population_size = POP_SIZE;
    settings.dimensions = DIM ;    
    settings.max_iterations = KMAX; 

       
    let lb =vec![-100.0f64; DIM];
    let ub =vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();  

    let mut fo = Sphere{};

    let mut algo : QAGO<Sphere> = QAGO::new(&settings, &mut fo);
    
    let result = algo.run();

    /*match result.convergence_trend{
        None => println!("QAGO: no convergence trend !!!"),
        Some(cv) => println!("QAGO: Convergence trend :\n {:?}", cv),
    };
    */

    /*match result.best_genome {
        None => println!("QAGO: no best solution !"),
        Some(bg)=> println!("QAGO: best-genome {:?}", bg),
    };
    */

    println!("QAGO : F1 (Sphere) test; Result: {:?}", result.to_string());
}

#[allow(dead_code)]
fn qago_f2_test1(){
    let mut settings : QAGOparams = QAGOparams::default();
    
    settings.population_size = POP_SIZE;
    settings.dimensions = DIM ;    
    settings.max_iterations = KMAX; 
       
    let lb =vec![-100.0f64; DIM];
    let ub =vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();  

    let mut fo = F2{};

    let mut algo : QAGO<F2> = QAGO::new(&settings, &mut fo);
    
    let result = algo.run();

    /*match result.convergence_trend{
        None => println!("QAGO: no convergence trend !!!"),
        Some(cv) => println!("QAGO: Convergence trend :\n {:?}", cv),
    };
    */

    /*match result.best_genome {
        None => println!("QAGO: no best solution !"),
        Some(bg)=> println!("QAGO: best-genome {:?}", bg),
    };
    */

    println!("QAGO : F2 test; Result: {:?}", result.to_string());
}

#[allow(dead_code)]
fn eo_f1_test1(){

    let mut settings : EOparams = EOparams::default();
    
    settings.population_size = POP_SIZE;
    settings.dimensions = DIM ;    
    settings.max_iterations = KMAX; 
    
    let lb =vec![-100.0f64; DIM];
    let ub =vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();    

    let mut fo = Sphere{};

    let mut eo : EO<Sphere> = EO::new(&settings, &mut fo);
    
    let result = eo.run();
       
    println!("EO result : \n best fitness : {:?} \n best genome : {:?}", result.best_fitness, result.best_genome);
    println!("Computation time : {:?}", result.computation_time);
    println!("Err : {:?}", result.err_report);
}


#[allow(dead_code)]
fn peo_f1_test1(){
    
    let mut settings : PSOparams = PSOparams::default();
    
    settings.population_size = POP_SIZE;
    settings.dimensions = DIM;    
    settings.max_iterations = KMAX; 
    
    let lb =vec![-100.0f64; DIM];
    let ub =vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();    

    let mut fo = Sphere{};

    let mut eo : PSO<Sphere> = PSO::new(&settings, &mut fo);
    
    let result = eo.run();
       
    println!("PSO result : \n best fitness : {:?} \n best genome : {:?}", result.best_fitness, result.best_genome);
    println!("Computation time : {:?}", result.computation_time);
    println!("Err : {:?}", result.err_report);
}

#[allow(dead_code)]
fn do_regression(){

    println!("------> REGRESSION PROBLEM :");
     
    let root = String::from("/home/sd/Documents/Rust_apps/sefar/src/bin/data");
      
    let path = format!("{}/{}", root, "Coxs_data_ALL.csv"); 

    //let mut settings : EOparams = EOparams::default();
    let mut settings : PSOparams = PSOparams::default();

    let dim = 14;

    settings.population_size = 80;
    settings.dimensions = dim;    
    settings.max_iterations = 2000; 
    
    let mut lb =vec![0.0f64; dim];
    lb[13] = 0.0f64;

    let ub =vec![10.0f64; dim];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();
    
    let mut regressn = Regression::new(path);

    //let mut eoa : EO<Regression> = EO::new(&settings, &mut regressn); 
    let mut eoa : PSO<Regression> = PSO::new(&settings, &mut regressn); 
   
    let result = eoa.run();
       
    println!("EO result : \n best fitness : {:?} \n best genome : {:?}", result.best_fitness, result.best_genome);
    println!("Computation time : {:?}", result.computation_time);
    println!("Err : {:?}", result.err_report);
    //println!("{:?}", result.convergence_trend);

    match result.best_genome {
        None => println!("error in optimization step..."),
        Some(genome) =>{

            let (rmsel, rmset, r2l, r2t) = regressn.compute_result_indexes(&genome.genes);

            println!("Indexes = {}, {}, {}, {}", rmsel, rmset, r2l, r2t);   
        },
    }

    
}

#[allow(dead_code)]
fn meo_test1(){

    let mut settings : EOparams = EOparams::default();
    
    settings.population_size = 7;
    settings.dimensions = 3 ;    
    settings.max_iterations = 2; 
    
    let lb =vec![-100.0f64; 3];
    let ub =vec![100.0f64; 3];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();    

    let mut fo = Sphere{};

    let mut eo : MEO<Sphere> = MEO::new(&settings, &mut fo);
    
    let result = eo.run();
       
    println!("EO result : \n best fitness : {:?} \n best genome : {:?}", result.best_fitness, result.best_genome);
    println!("Computation time : {:?}", result.computation_time);
    println!("Err : {:?}", result.err_report);

}

#[allow(dead_code)]
fn para_meo_test1(){

    let mut settings : EOparams = EOparams::default();
    
    settings.population_size = 7;
    settings.dimensions = 3 ;    
    settings.max_iterations = 2; 
    
    let lb =vec![-100.0f64; 3];
    let ub =vec![100.0f64; 3];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();    

    let mut fo = Sphere{};

    let mut eo : MEO<Sphere> = MEO::new(&settings, &mut fo);
    
    let result = eo.run();
       
    println!("EO result : \n best fitness : {:?} \n best genome : {:?}", result.best_fitness, result.best_genome);
    println!("Computation time : {:?}", result.computation_time);
    println!("Err : {:?}", result.err_report);

}

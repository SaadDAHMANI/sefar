use sefar::core::eoa::EOA;
use sefar::benchmarks::functions::{Sphere, F2};
use sefar::core::optimization_result::OptimizationResult;
use sefar::algos::eo::{EO, EOparams};
use sefar::algos::pso::{PSO, PSOparams};

use sefar::algos::meo::MEO;
// use sefar::algos::qago::{QAGOparams, QAGO};
use sefar::algos::go::{GOparams, GO};

const DIM : usize = 10;
const POP_SIZE : usize = 30;
const KMAX : usize = 200; 

fn main() {
    println!("Hello, sefar !");
   
    println!("Evaluation with Max_Iter = {}", KMAX);

    println!("______________________GO : F1______________________");

    //go_f1_test1();
   
    /*
    println!("______________________GO : F2______________________");
    go_f2_test1();
   */
    //println!("_______________________QAGO : F1______________________");

    //qago_f1_test1();

    //pco_f1_test1();
    
    //eo_f1_test1();

    //peo_f1_test1();

    //meo_test1();

    #[cfg(feature ="binary")] {
        println!("Run Binary tests");
        eo_f1_binary_test();
        println!("_________________________________________________________________");                                 
        go_f1_binary_test();
    }
    
  }

  ///
  /// run the binary version of Growth Optimizer (Binary-GO).
  /// 
  #[cfg(feature = "binary")]
  #[allow(dead_code)]
  fn go_f1_binary_test(){

    // Define the parameters of GO:
    let search_agents : usize = POP_SIZE;
    let dim : usize = DIM;
    let max_iterations : usize = KMAX;
    let lb = vec![0.0; DIM];
    let ub = vec![1.0; DIM];
    
    // Build the parameter struct:
    let settings : GOparams = GOparams::new(search_agents, dim, max_iterations, &lb, &ub);
    
    // Define the problem to optimize:
    let mut fo = Sphere{};
  
    // Build the optimizer:
    let mut algo : GO<Sphere> = GO::new(&settings, &mut fo);
    
    // Run the GO algorithm: 
    let result : OptimizationResult = algo.run();

    // Print the results:
    println!("The optimization results of Binary-GO : {}", result.to_string());

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
  
/* #[allow(dead_code)]
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
} */

/* #[allow(dead_code)]
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
 */
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

#[cfg(feature="binary")]
#[allow(dead_code)]
fn eo_f1_binary_test(){

    let mut settings : EOparams = EOparams::default();
    
    settings.population_size = POP_SIZE;
    settings.dimensions = DIM;    
    settings.max_iterations = KMAX; 
    
    let lb =vec![0.0f64; DIM];
    let ub =vec![1.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();    

    let mut fo = Sphere{};

    let mut eo : EO<Sphere> = EO::new(&settings, &mut fo);
    
    let result = eo.run();
       
    println!("Binary-EO result : \n best fitness : {:?} \n best genome : {:?}", result.best_fitness, result.best_genome);
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

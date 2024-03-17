include!("strm_regression.rs");
include!("dataset.rs");

use sefar::core::eoa::EOA;
use sefar::benchmarks::functions::Sphere;
use sefar::sequential_algos::eo::{EO, EOparams};
use sefar::sequential_algos::pso::{PSO, PSOparams};
use sefar::sequential_algos::pco::{PCO, PCOparams};

use sefar::sequential_algos::meo::MEO;
use sefar::sequential_algos::qago::{QAGOparams, QAGO};

const DIM : usize = 3;
const POP_SIZE : usize = 4;
const KMAX : usize = 1;



fn main() {
    println!("Hello, sefar !");

    //qago_f1_test1();

    //pco_f1_test1();
    
    //eo_f1_test1();

    //println!("_______________________________________________");

    //peo_f1_test1();

    //meo_test1();

    //do_regression();
    
    let fitness = [1.0, 5.0, 3.0, 4.0];
                                             
    let mut ind : Vec<usize> = (0..fitness.len()).collect();
    ind.sort_by(|&a, &b| fitness[a].partial_cmp(&fitness[b]).unwrap());

    println!("indexes : {:?}", ind);     


     let x = [10, 5, 3, 1];
     let z = match x.iter().enumerate().min_by_key(|&(_, v)| v) {
        Some((i, _))=> i,
        None =>0,   
     };

       println!("The best index is  z = {}", z);                                                                                         

    let x : f64 = -2.7; 
    let y= x.ceil();
    println!("ceil x : {}", y);                                            
                                             

  }



#[allow(dead_code)]
fn qago_f1_test1(){
    let mut settings : QAGOparams = QAGOparams::default();
    
    /* settings.population_size = POP_SIZE;
    settings.dimensions = DIM ;    
    settings.max_iterations = KMAX; 

       
    let lb =vec![-100.0f64; DIM];
    let ub =vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();  
 */
    let mut fo = Sphere{};

    let mut algo : QAGO<Sphere> = QAGO::new(&settings, &mut fo);
    
    let result = algo.run();

    match result.convergence_trend{
        None => println!("QAGO: no convergence trend !!!"),
        Some(cv) => println!("QAGO: Convergence trend :\n {:?}", cv),
    };

}

#[allow(dead_code)]
fn pco_f1_test1(){
    let mut settings : PCOparams = PCOparams::default();
    
    settings.population_size = POP_SIZE;
    settings.dimensions = DIM ;    
    settings.max_iterations = KMAX; 

       
    let lb =vec![-100.0f64; DIM];
    let ub =vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();  

    let mut fo = Sphere{};

    let mut pco : PCO<Sphere> = PCO::new(&settings, &mut fo);
    
    let result = pco.run();

    match result.convergence_trend{
        None => println!("PCO: no convergence trend !!!"),
        Some(cv) => println!("PCO: Convergence trend :\n {:?}", cv),
    };

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

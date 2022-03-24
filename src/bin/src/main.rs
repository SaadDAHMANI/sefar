
use sefar::benchmarks::functions::F1;
use sefar::sequential_algos::eo::EO;
use sefar::sequential_algos::eo::EOparams;
use sefar::core::eoa::EOA;

fn main() {
    println!("Hello, sefar !");
     
    fn f1(x : i32)-> i32 {
        x+10i32
    }

    let a : i32 = 23;
    let b = f1(a);
    println!("b : {:?}", b);

    //eo_f1_test1();
}

#[allow(dead_code)]
fn eo_f1_test1(){
    
    let dim : usize =30;

    let mut settings : EOparams = EOparams::default();
    
    settings.population_size = 30;
    settings.dimensions = dim;    
    settings.max_iterations = 500; 
    
    let lb =vec![-100.0f64; dim];
    let ub =vec![100.0f64; dim];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();    

    let mut fo = F1 {};

    let mut eo : EO<F1> = EO::new(&settings, &mut fo);
    
    let result = eo.run();
       
    println!("EO result : \n best fitness : {:?} \n best genome : {:?}", result.best_fitness, result.best_genome);
    println!("Computation time : {:?}", result.computation_time);
    println!("Err : {:?}", result.err_report);
}
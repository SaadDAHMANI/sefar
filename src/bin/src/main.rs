
use sefar::benchmarks::functions::Sphere;
use sefar::sequential_algos::eo::EO;
use sefar::sequential_algos::eo::EOparams;
use sefar::sequential_algos::pso::PSO;
use sefar::sequential_algos::pso::PSOparams;
use sefar::core::eoa::EOA;

const DIM : usize =10;
const POP_SIZE : usize =30;
const KMAX : usize = 500;

fn main() {
    println!("Hello, sefar !");
    
    eo_f1_test1();

    println!("_______________________________________________");

    peo_f1_test1();
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




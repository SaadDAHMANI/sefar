
use sefar::benchmarks::functions::F1;
use sefar::sequential_algos::eo::EO;
use sefar::sequential_algos::eo::EOparams;
use sefar::core::eoa::EOA;

fn main() {
    println!("Hello, sefar !");

    eo_f1_test1();
}

fn eo_f1_test1(){
    let settings : EOparams = EOparams::default();
    //settings.population_size = 20;
    //settings.max_iterations = 500; 

    let mut fo = F1 {};

    let mut eo : EO<F1> = EO::new(&settings, &mut fo);
    
    let result = eo.run();
       
    println!("EO result : \n best fitness : {:?} \n best genome : {:?}", result.best_fitness, result.best_genome);

 }
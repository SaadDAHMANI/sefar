use sefar::algos::eo::{BiEO, BiEOparams, EOparams, EO, MEO};

use sefar::algos::pso::{PSOparams, PSO};
use sefar::benchmarks::functions::{Sphere, SumAbsFunction, F2};
use sefar::core::eoa::EOA;
use sefar::core::optimization_result::OptimizationResult;

// use sefar::algos::qago::{QAGOparams, QAGO};
// use sefar::algos::apgsk::{APGSKparams, APGSK};

use sefar::algos::go::{GOparams, GO};
use sefar::algos::gsk::{GSKparams, GSK};

//use sefar::algos::lshade_spacma::{LshadeSpacma, LshadeSpacmaParams};

const DIM: usize = 5;
const POP_SIZE: usize = 12;
const KMAX: usize = 10;

fn main() {
    // lshade_spacma_test1();

    //apgsk_f1_test1();

    // #[cfg(not(feature = "parallel"))]
    gsk_f1_test1();

    //--------------------------------------------------------------------
    //println!("Evaluation with Max_Iter = {}", KMAX);

    // println!("______________________GO : F1______________________");
    // go_f1_test1();

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
    //bieo_f1_binary_test();

    println!("------------------GSK : F1---------------------");
    gsk_f1_test1();
}
/*
#[allow(dead_code)]
fn lshade_spacma_test1() {
    let mut settings: LshadeSpacmaParams = LshadeSpacmaParams::default();
    settings.population_size = POP_SIZE;
    settings.problem_dimension = DIM;
    settings.max_iterations = KMAX;
    let lb = vec![-10.0f64; DIM];
    let ub = vec![10.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();

    let mut fo: SumAbsFunction = SumAbsFunction {};
    let mut algo: LshadeSpacma<SumAbsFunction> = LshadeSpacma {
        problem: &mut fo,
        params: &settings,
    };

    let result = algo.run();
    println!(
        "LSHADE_SPACMA : F0 (SumAbsFunction) test; Result: {:?}",
        result.to_string()
    );
}
*/

#[allow(dead_code)]
fn apgsk_f1_test1() {
    let mut settings: GSKparams = GSKparams::default(); // APGSKparams = APGSKparams::default();

    settings.population_size = POP_SIZE;
    settings.problem_dimension = DIM;
    settings.max_iterations = KMAX;

    let lb = vec![-100.0f64; DIM];
    let ub = vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();

    settings.partition_size_p = 0.10;

    let mut fo = SumAbsFunction {}; // Sphere{};

    let mut algo: GSK<SumAbsFunction> = GSK::new(&settings, &mut fo); //APGSK<SumAbsFunction> = APGSK::new(&settings, &mut fo);

    let result: OptimizationResult = algo.run();

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

    println!(
        "Gaining-Sharing Knowledge optimizer with Adaptive Parameters (APGSK) : F0 (SumAbsFunction) test; Result: {:?}",
        result.to_string()
    );
}

#[allow(dead_code)]
fn gsk_f1_test1() {
    let mut settings: GSKparams = GSKparams::default();

    settings.population_size = POP_SIZE;
    settings.problem_dimension = DIM;
    settings.max_iterations = KMAX;

    let lb = vec![-100.0f64; DIM];
    let ub = vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();

    settings.partition_size_p = 0.2;
    settings.kr = 0.7;

    let mut fo = Sphere {};

    let mut algo: GSK<Sphere> = GSK::new(&settings, &mut fo);

    let result: OptimizationResult = algo.run();

    println!(
        "Gaining-Sharing Knowledge optimizer (GSK) : F1 (Sphere) test; Result: {:?}",
        result.to_string()
    );
}

#[cfg(feature = "parallel")]
#[allow(dead_code)]
fn para_gsk_f1_test1() {
    let mut settings: GSKparams = GSKparams::default();

    settings.population_size = POP_SIZE;
    settings.problem_dimension = DIM;
    settings.max_iterations = KMAX;

    let lb = vec![-100.0f64; DIM];
    let ub = vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();

    //settings.partition_size_p = 0.2;
    settings.kr = 0.8;

    let mut fo = Sphere {};

    let mut algo: GSK<Sphere> = GSK::new(&settings, &mut fo);

    let result: OptimizationResult = algo.run();

    println!(
        "Parallel Gaining-Sharing Knowledge optimizer (Para-GSK) : F1 (Sphere) test; Result: {:?}",
        result.to_string()
    );
}

#[allow(dead_code)]
fn go_f1_test1() {
    let mut settings: GOparams = GOparams::default();

    settings.population_size = POP_SIZE;
    settings.problem_dimension = DIM;
    settings.max_iterations = KMAX;

    let lb = vec![-100.0f64; DIM];
    let ub = vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();

    let mut fo = Sphere {};

    let mut algo: GO<Sphere> = GO::new(&settings, &mut fo);

    let result: OptimizationResult = algo.run();
    println!("The optimization results of GO : {}", result.to_string());

    // Save the results:
    let file = "/home/sd/Documents/Rust_apps/GO.csv";

    let header = format!("{}", settings);

    let _result = result.save(Some(&header), &file);
    println!("Saving optimization results in teh file: {:?}", _result);
}

#[allow(dead_code)]
fn go_f2_test1() {
    let mut settings: GOparams = GOparams::default();

    settings.population_size = POP_SIZE;
    settings.problem_dimension = DIM;
    settings.max_iterations = KMAX;

    let lb = vec![-100.0f64; DIM];
    let ub = vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();

    let mut fo = F2 {};

    let mut algo: GO<F2> = GO::new(&settings, &mut fo);

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

    println!(
        "Growth Optimizer (GO) : F2 test; Result: {:?}",
        result.to_string()
    );
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
fn eo_f1_test1() {
    let mut settings: EOparams = EOparams::default();

    settings.population_size = POP_SIZE;
    settings.problem_dimension = DIM;
    settings.max_iterations = KMAX;

    let lb = vec![-100.0f64; DIM];
    let ub = vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();

    let mut fo = Sphere {};

    let mut eo: EO<Sphere> = EO::new(&settings, &mut fo);

    let result = eo.run();

    // Print the results:
    println!("The optimization results of GO : {}", result.to_string());
}

#[allow(dead_code)]
fn bieo_f1_binary_test() {
    let mut settings: BiEOparams = BiEOparams::default();

    settings.population_size = POP_SIZE;
    settings.problem_dimension = DIM;
    settings.max_iterations = KMAX;

    let mut fo = Sphere {};

    let mut bieo: BiEO<Sphere> = BiEO::new(&settings, &mut fo);

    let result = bieo.run();

    println!(
        "Binary-EO result : \n best fitness : {:?} \n best genome : {:?}",
        result.best_fitness, result.best_genome
    );
    println!("Computation time : {:?}", result.computation_time);
    println!("Err : {:?}", result.err_report);
}

#[allow(dead_code)]
fn peo_f1_test1() {
    let mut settings: PSOparams = PSOparams::default();

    settings.population_size = POP_SIZE;
    settings.problem_dimension = DIM;
    settings.max_iterations = KMAX;

    let lb = vec![-100.0f64; DIM];
    let ub = vec![100.0f64; DIM];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();

    let mut fo = Sphere {};

    let mut eo: PSO<Sphere> = PSO::new(&settings, &mut fo);

    let result = eo.run();

    println!(
        "PSO result : \n best fitness : {:?} \n best genome : {:?}",
        result.best_fitness, result.best_genome
    );
    println!("Computation time : {:?}", result.computation_time);
    println!("Err : {:?}", result.err_report);
}

#[allow(dead_code)]
fn meo_test1() {
    let mut settings: EOparams = EOparams::default();

    settings.population_size = 7;
    settings.problem_dimension = 3;
    settings.max_iterations = 2;

    let lb = vec![-100.0f64; 3];
    let ub = vec![100.0f64; 3];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();

    let mut fo = Sphere {};

    let mut eo: MEO<Sphere> = MEO::new(&settings, &mut fo);

    let result = eo.run();

    println!(
        "EO result : \n best fitness : {:?} \n best genome : {:?}",
        result.best_fitness, result.best_genome
    );
    println!("Computation time : {:?}", result.computation_time);
    println!("Err : {:?}", result.err_report);
}

#[allow(dead_code)]
fn para_meo_test1() {
    let mut settings: EOparams = EOparams::default();

    settings.population_size = 7;
    settings.problem_dimension = 3;
    settings.max_iterations = 2;

    let lb = vec![-100.0f64; 3];
    let ub = vec![100.0f64; 3];

    settings.lower_bounds = lb.as_slice();
    settings.upper_bounds = ub.as_slice();

    let mut fo = Sphere {};

    let mut eo: MEO<Sphere> = MEO::new(&settings, &mut fo);

    let result = eo.run();

    println!(
        "EO result : \n best fitness : {:?} \n best genome : {:?}",
        result.best_fitness, result.best_genome
    );
    println!("Computation time : {:?}", result.computation_time);
    println!("Err : {:?}", result.err_report);
}

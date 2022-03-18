pub mod parameters;
pub mod problem;
pub mod genome;
pub mod optimization_result;
pub mod eoa;

use optimization_result::OptimizationResult;
pub trait EvolutionaryAlgo {
    fn run_epoch(&mut self){}
    fn run(&mut self)-> Option<OptimizationResult> {
        None
    }
}
 
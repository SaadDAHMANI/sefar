use crate::core::optimization_result::OptimizationResult;

///
/// Public trait for Evolutionary Optimization Algorithms
/// 
pub trait EOA {
    
    ///
    /// Run algorithm until reach stopping criterion
    /// 
    fn run(&mut self)-> OptimizationResult {
        let  result = OptimizationResult{
             best_genome : None,
             best_fitness :None,
             convergence_trend : None,
             computation_time : None,
             err_report : None, 
        };
        result 
    }
    
}
//
// The Objectivefunction can be used with sequential algorithms (mutation of problem can be done).
//
use crate::core::genome::Genome;
pub trait Problem: Send + Sync + Clone {
    ///
    /// Define the objective function to be called in sequential mode.
    /// In sequential run, the problem can be modified (safely) when the objective function is called.
    ///
    #[cfg(not(feature = "parallel"))]
    fn objectivefunction(&mut self, genome: &mut [f64]) -> f64 {
        genome.iter().fold(0.0f64, |sum, x| sum + x)
    }

    ///
    /// Define a custom behavior when ierations progress.
    ///
    #[allow(dead_code)]
    fn iteration_increment(
        &self,
        _current_iteration: usize,
        _current_best_genome: &Genome,
        _break_process: &mut bool,
    ) {
    }
    ///
    /// Define the objective function to be called in parallel mode.
    /// In parallel mode, the problem cannot be modified when the objective function is called.
    /// In this version, multiple threads cannot safely modify the struct that implements the problem trait.
    ///
    #[cfg(feature = "parallel")]
    fn objectivefunction(&self, genome: &[f64]) -> f64 {
        genome.iter().fold(0.0f64, |sum, x| sum + x)
    }
}

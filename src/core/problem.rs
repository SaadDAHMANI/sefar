//
// The Objectivefunction can be used with sequential algorithms (mutation of problem can be done).
//
use crate::core::genome::Genome;
pub trait Problem: Send + Sync + Clone {
    ///
    /// Define the objective function to be called in sequential mode.
    /// In sequential run, the problem can be modified (safely) when the objective function is called.
    ///
    /// # Arguments
    ///
    /// * `gemone` - The candidate solution. Use it to compute it fitness.
    ///
    /// # Returns
    ///
    /// * The fitness value of the `genome` as f64.
    ///  
    #[cfg(not(feature = "parallel"))]
    fn objectivefunction(&mut self, genome: &mut [f64]) -> f64 {
        genome.iter().fold(0.0f64, |sum, x| sum + x.abs())
    }

    ///
    /// Define a custom behavior when ierations progress.
    ///
    /// # Arguments
    ///
    /// * `current_iter` - The current iteration.
    /// * `current_best_genome` - The current best solution.
    /// * `beak_process` - Break optimization process by :`*break_process = true`;
    ///
    #[allow(dead_code)]
    #[allow(unused_variables)]
    fn iteration_increment(
        &self,
        current_iter: usize,
        current_best_genome: &Genome,
        break_process: &mut bool,
    ) {
    }
    ///
    /// Define the objective function to be called in parallel mode.
    /// In parallel mode, the problem cannot be modified when the objective function is called.
    /// In this version, multiple threads cannot safely modify the struct that implements the problem trait.
    ///
    /// # Arguments
    ///
    /// * `gemone` - The candidate solution. Use it to compute it fitness.
    ///
    /// # Returns
    ///
    /// * The fitness value of the `genome` as f64.
    #[cfg(feature = "parallel")]
    fn objectivefunction(&self, genome: &[f64]) -> f64 {
        genome.iter().fold(0.0f64, |sum, x| sum + x)
    }
}

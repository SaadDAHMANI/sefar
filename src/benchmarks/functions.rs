use crate::core::problem::Problem;
//use crate::core::genome::Genome;

#[derive(Debug, Clone)]
pub struct SumAbsFunction {}
impl Problem for SumAbsFunction {
    #[cfg(not(feature = "parallel"))]
    fn objectivefunction(&mut self, genome: &[f64]) -> f64 {
        let fitness = genome.iter().fold(0.0f64, |sum, g| sum + g.abs());
        fitness
    }

    fn iteration_increment(
        &self,
        current_iter: usize,
        current_best_genome: &crate::core::genome::Genome,
        break_process: &mut bool,
    ) {
        println!(
            "Iter: {}, current best-fit : {:?}",
            current_iter, current_best_genome.fitness
        );

        if current_iter > 100 {
            *break_process = true;
        }
    }

    #[cfg(feature = "parallel")]
    fn objectivefunction(&self, genome: &[f64]) -> f64 {
        let fitness = genome.iter().fold(0.0f64, |sum, g| sum + g.abs());
        fitness
    }
}

///
/// Sphere benchmark function (F1).
/// Fi(X) = Sum(|X|)
/// where X = {x1, x2, ..... xd}, and 'd' is the problem dimension.
///
#[derive(Debug, Clone)]
pub struct Sphere {}

impl Problem for Sphere {
    ///
    /// Define the objective function. The later is called in sequential mode.
    ///
    #[cfg(not(feature = "parallel"))]
    fn objectivefunction(&mut self, genome: &[f64]) -> f64 {
        let fitness = genome.iter().fold(0.0f64, |sum, g| sum + g.powi(2));
        //let fitness = genome.iter().fold(0.0f64, |sum, g| sum + g);
        fitness
    }

    #[cfg(feature = "parallel")]
    fn objectivefunction(&self, genome: &[f64]) -> f64 {
        let fitness = genome.iter().fold(0.0f64, |sum, g| sum + g.powi(2));
        fitness
    }

    fn iteration_increment(
        &self,
        current_iteration: usize,
        current_best_genome: &crate::core::genome::Genome,
        break_process: &mut bool,
    ) {
        if let Some(best_fit) = current_best_genome.fitness {
            if best_fit > 0.00001 {
                println!(
                    "Iteration : {}, best.fitness : {:?}",
                    current_iteration, best_fit
                );
            } else {
                *break_process = true;
            }
        }
    }
}

///
/// F2 = Sum(x_i) + Prod(x_i) = Sum(X) + Prod(X).
/// Search space : [-100.0, 100.0];
/// Problem dimension : 30;
/// Optimal value = vec![0.0; 30];
#[derive(Debug, Clone)]
pub struct F2 {}
impl Problem for F2 {
    #[cfg(not(feature = "parallel"))]
    fn objectivefunction(&mut self, genome: &[f64]) -> f64 {
        let sum = genome.iter().fold(0.0f64, |sum, g| sum + g.abs());
        let prod = genome.iter().fold(1.0f64, |prod, g| prod * f64::abs(*g));
        sum + prod
    }

    #[cfg(feature = "parallel")]
    fn objectivefunction(&self, genome: &[f64]) -> f64 {
        let sum = genome.iter().fold(0.0f64, |sum, g| sum + g.abs());
        let prod = genome.iter().fold(1.0f64, |prod, g| prod * f64::abs(*g));
        sum + prod
    }
}

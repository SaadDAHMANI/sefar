use crate::core::problem::Problem;
//use crate::core::genome::Genome;

#[derive(Debug,Clone)]

///
/// Sphere enchmark function
/// 
pub struct F1{}

impl Problem for F1{
    fn objectivefunction(&mut self, genome : &[f64]) ->f64 {
        let fitness = genome.iter().fold(0.0f64, |sum, g| sum + g.powi(2));
        fitness
    }
}

#[derive(Debug, Clone)]
pub struct F2{}
impl Problem for F2{
    fn objectivefunction(&mut self, genome : &[f64]) ->f64 {
      let sum = genome.iter().fold(0.0f64, |sum, g| sum + g.abs());
      let prod = genome.iter().fold(1.0f64, |prod, g| prod * g.abs());
      sum + prod
    }
}

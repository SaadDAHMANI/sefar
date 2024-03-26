use crate::core::problem::Problem;
//use crate::core::genome::Genome;

#[derive(Debug,Clone)]

///
/// Sphere enchmark function (F1)
/// 
pub struct Sphere{}

impl Problem for Sphere{
    fn objectivefunction(&mut self, genome : &[f64]) ->f64 {
        let fitness = genome.iter().fold(0.0f64, |sum, g| sum + g.powi(2));
        fitness        
    }
}

///
/// F2 = Sum(x_i) + Prod(x_i) = Sum(X) + Prod(X)  
/// 
#[derive(Debug, Clone)]
pub struct F2{}
impl Problem for F2 {
    fn objectivefunction(&mut self, genome : &[f64]) ->f64 {
      let sum = genome.iter().fold(0.0f64, |sum, g| sum + g.abs());
      let prod = genome.iter().fold(1.0f64, |prod, g| prod* f64::abs(*g));
      sum + prod
    }
}

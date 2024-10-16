use crate::core::problem::Problem;
//use crate::core::genome::Genome;


///
/// Sphere benchmark function (F1). 
/// Fi(X) = Sum(|X|)
/// where X = {x1, x2, ..... xd}, and 'd' is the problem dimension.
/// 
#[derive(Debug,Clone)]
pub struct Sphere{}

impl Problem for Sphere{
  ///
  /// Define the objective function. The later is called in sequential mode. 
  ///   
  #[cfg(not(feature="parallel"))]
  fn objectivefunction(&mut self, genome : &[f64]) ->f64 {
    let fitness = genome.iter().fold(0.0f64, |sum, g| sum + g.powi(2));
      fitness        
  }
  
  #[cfg(feature="parallel")]
  fn objectivefunction(&self, genome : &[f64]) ->f64 {
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
    #[cfg(not(feature="parallel"))]
    fn objectivefunction(&mut self, genome : &[f64]) ->f64 {
      let sum = genome.iter().fold(0.0f64, |sum, g| sum + g.abs());
      let prod = genome.iter().fold(1.0f64, |prod, g| prod* f64::abs(*g));
      sum + prod
    }

    #[cfg(feature="parallel")]
    fn objectivefunction(&self, genome : &[f64]) ->f64 {
      let sum = genome.iter().fold(0.0f64, |sum, g| sum + g.abs());
      let prod = genome.iter().fold(1.0f64, |prod, g| prod* f64::abs(*g));
      sum + prod
    }
} 

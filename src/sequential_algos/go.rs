use crate::core::eoa::EOA;
use crate::core::genome::Genome;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;




#[derive(Debug, Clone)]
pub struct GOparams<'a> {
    pub population_size : usize,
    pub dimensions : usize,
    pub max_iterations: usize,
    pub lower_bounds : &'a [f64],
    pub upper_bounds : &'a [f64],
}

impl<'a> Parameters for GOparams<'a> {
    fn get_dimensions(&self)->usize {
         self.dimensions    
    }

    fn get_max_iterations(&self)->usize {
        usize::max(self.max_iterations,1) 
    }

    fn get_population_size(&self)->usize {
        usize::max(self.population_size, 5)
    }

    fn get_lower_bounds(&self)->Vec<f64> {
        self.lower_bounds.to_vec()
    }

    fn get_upper_bounds(&self)->Vec<f64> {
        self.upper_bounds.to_vec()
    }
}

impl<'a> Default for GOparams<'a>{

    ///
    /// Return default values of parameters, as following :
    /// 
    /// ~~~
    /// 
    ///  use sefar::sequential_algos::go::*;
    /// 
    ///  GOparams{
    ///     population_size : 10,
    ///     dimensions : 3,
    ///     max_iterations : 100,
    ///     lower_bounds : &[100.0f64, 100.0, 100.0],
    ///     upper_bounds : &[-100.0f64, -100.0, -100.0],
    /// };
    /// ~~~
    /// 
    fn default()->Self{
        GOparams {
            population_size : 10,
            dimensions : 3,
            max_iterations : 1,
            lower_bounds : &[-100.0f64, -100.0, -100.0],
            upper_bounds : &[100.0f64, 100.0, 100.0],
        }
    }
}


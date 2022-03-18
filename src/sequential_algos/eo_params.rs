use crate::core::parameters::Parameters;

#[derive(Debug, Clone)]
pub struct EOparams<'a> {
    pub population_size : usize,
    pub dimensions : usize,
    pub max_iterations: usize,
    pub lower_bounds : &'a [f64],
    pub upper_bounds : &'a [f64],
    pub a1 : f64,
    pub a2 : f64,
    pub gp : f64,
}

#[allow(dead_code)]
impl<'a> EOparams<'a>{
    pub fn new(p_size: usize, dim : usize, max_iter : usize, lb : &'a [f64], 
    ub : &'a [f64], a1 :f64, a2 :f64, gp : f64)->EOparams<'a> {
        EOparams{
            population_size : p_size,
            dimensions : dim,
            max_iterations : max_iter,
            lower_bounds : lb,
            upper_bounds : ub,
            a1,
            a2,
            gp,
        }
    }    
}

impl<'a> Parameters for EOparams<'a> {

    fn get_population_size(&self)->usize{
        self.population_size
    }

    fn get_dimensions(&self)-> usize {
        self.dimensions
    }

    fn get_max_iterations(&self)->usize{
        self.max_iterations
    }

    fn get_lower_bounds(&self)-> Vec<f64>{
        self.lower_bounds.to_vec()
    }

    fn get_upper_bounds(&self)-> Vec<f64>{
        self.upper_bounds.to_vec()
    }        
}





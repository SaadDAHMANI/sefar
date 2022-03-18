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
    
    ///
    /// Return default values of parameters, as following :
    /// 
    /// ~~~
    /// 
    ///  use sefar::sequential_algos::eo_params::EOparams;
    /// 
    ///  EOparams{
    ///     population_size : 10,
    ///     dimensions : 3,
    ///     max_iterations : 100,
    ///     lower_bounds : &[100.0f64, 100.0, 100.0],
    ///     upper_bounds : &[-100.0f64, -100.0, -100.0],
    ///     a1 : 2.0f64,
    ///     a2 : 1.0f64,
    ///     gp : 0.5f64,
    /// };
    /// ~~~
    /// 
    pub fn default()->Self{
        EOparams{
            population_size : 10,
            dimensions : 3,
            max_iterations : 100,
            lower_bounds : &[100.0f64, 100.0, 100.0],
            upper_bounds : &[-100.0f64, -100.0, -100.0],
            a1 : 2.0f64,
            a2 : 1.0f64,
            gp : 0.5f64,
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

#[cfg(test)]
mod eo_params_tests {
    
    use super::*;

    #[test]
    fn test_ub_slice(){
        let d : usize = 5;
        let n : usize =10;
        let k : usize = 100;
        
        let ub = vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let lb = ub.clone();
        
        let params = EOparams{
            population_size : n,
            max_iterations : k,
            dimensions : d,
            lower_bounds : lb.as_slice(),
            upper_bounds : ub.as_slice(),
            a1 : 2.0f64,
            a2 : 1.0f64,
            gp : 0.5f64,
        };

        let sl_ub =  vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let slice_ub = sl_ub.as_slice();
        
        assert_eq!(params.upper_bounds, slice_ub);
    }

    #[test]
    fn test_default_fn(){
        let p = EOparams::default();
        assert_eq!(p.a1, 2.0f64);
        assert_eq!(p.a2, 1.0f64);
        assert_eq!(p.gp, 0.50f64);
    }
}



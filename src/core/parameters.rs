
pub trait Parameters {
   
    fn get_population_size(&self)->usize {
          10usize
    }
    
    fn get_dimensions(&self)->usize {
        10usize
    }

    fn get_max_iterations(&self)->usize{
        1usize
    }

    fn get_lower_bounds(&self)->Vec<f64>{
        let d = self.get_dimensions();
        let mut lb = Vec::new();
        for _i in 0..d{
            lb.push(-100.0f64);
        } 
        lb
    }

    fn get_upper_bounds(&self)->Vec<f64>{
        let d = self.get_dimensions();
        let mut ub = Vec::new();
        for _i in 0..d{
            ub.push(100.0f64);
        } 
        ub
    }   
}  
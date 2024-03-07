//
// Implementation of Equilibrium Optimizer (EO)
// 
 
extern crate rand;
use rand::distributions::{Distribution, Uniform};
//use rand::prelude::ThreadRng;
use std::time::Instant;

use crate::core::eoa::EOA;
use crate::core::genome::Genome;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
use crate::core::optimization_result::OptimizationResult;
use crate::common::*;

///
/// Sequential Equilibrium Optimizer (EO)
/// Reference:
/// "Faramarzi, A., Heidarinejad, M., Stephens, B., & Mirjalili, S. (2020).
/// Equilibrium optimizer: A novel optimization algorithm. Knowledge-Based Systems, 191, 105190."
/// 
#[derive(Debug)]
pub struct EO<'a, T : Problem> {
     pub problem : &'a mut T,
     pub params : &'a EOparams<'a>,
     pub optimization_result : OptimizationResult,
}

impl<'a, T : Problem> EO<'a, T>{

    pub fn new(settings :&'a EOparams, problem : &'a mut T )->Self{       
        let result = OptimizationResult{
            best_genome : None,
            best_fitness :None,
            convergence_trend : None,
            computation_time : None,
            err_report : None, 
        };
       
        EO { 
             problem,
             params: settings,
             optimization_result: result,            
        }
    }   
}

impl<'a, T: Problem> EOA for EO<'a, T> {
   
    fn run(&mut self)-> OptimizationResult{
        
        let chronos = Instant::now();
    
        //check paramaters
        //let params = self.params.clone();
        
        match self.params.check(){
            Err(error) => OptimizationResult::get_none(error),
            Ok(()) =>{
                let dim =  self.params.get_dimensions();
                let particles_no = self.params.get_population_size();
                let lb = self.params.get_lower_bounds();
                let ub = self.params.get_upper_bounds();
                let max_iter = self.params.get_max_iterations();
    
                //
                // a1=2;
                // a2=1;
                // GP=0.5;
                let a1 : f64 = self.params.a1;
                let a2 : f64 = self.params.a2;
                let gp : f64 = self.params.gp;
    
                // Initialize variables 
                //Ceq1=zeros(1,dim);   Ceq1_fit=inf; 
                //Ceq2=zeros(1,dim);   Ceq2_fit=inf; 
                //Ceq3=zeros(1,dim);   Ceq3_fit=inf; 
                //Ceq4=zeros(1,dim);   Ceq4_fit=inf;
                
                let mut ceq1 = vec![0.0f64; dim];
                
                let mut ceq2 = vec![0.0f64; dim];
                
                let mut ceq3 = vec![0.0f64; dim];
                
                let mut ceq4 = vec![0.0f64; dim];
                
                let mut ceq_ave = vec![0.0f64; dim];
                
                let mut ceq1_fit = f64::MAX;
                
                let mut ceq2_fit = f64::MAX;
                
                let mut ceq3_fit = f64::MAX;
                
                let mut ceq4_fit = f64::MAX;    
                
                let mut ceq1_index : usize = 0;
                //let mut ceq2_index : usize = 0;
                //let mut ceq3_index : usize = 0;
                //let mut ceq4_index : usize = 0;
        
                // Iter=0; V=1;
                let mut iter =0;
                let v : f64 = 1.0;
                        
                // to store agents fitness values
                let mut fitness = vec![0.0f64; particles_no];
                let mut fit_old = vec![0.0f64; particles_no];
                let mut c_old = vec![vec![0.0f64; dim]; particles_no];
                let mut c_pool = vec![vec![0.0f64; dim]; 5];
    
                let mut lambda = vec![0.0f64; dim];
                let mut r = vec![0.0f64; dim];
                let mut r1 = vec![0.0f64; dim];
                let mut r2 = vec![0.0f64; dim];
                let mut ceq = vec![0.0f64; dim];
                let mut f = vec![0.0f64; dim];
                let mut _gcp :f64 =0.0;
                //------------------------------------------
                let interval = Uniform::from(0..c_pool.len());
                //let between01 = Uniform::from(0.0..=1.0);
                let mut rng = rand::thread_rng();
                //------------------------------------------
                
                let mut convergence_curve = vec![0.0f64; max_iter]; 
                let mut _index : usize = 0;
                let mut _g0 : f64 = 0.0; 
                let mut _g : f64 = 0.0;
                
                //let chronos = Instant::now();
                
                //C=initialization(Particles_no,dim,ub,lb);
                let mut c = self.initialize(self.params);
        
                // the main loop of EO
                while iter < max_iter {
                
                    for i in 0..c.len() {
                
                        // space bound
                        for j in 0..dim {
                            if  c[i].genes[j] < lb[j] { c[i].genes[j] = lb[j];}
                            
                            if c[i].genes[j] > ub[j] { c[i].genes[j] = ub[j];}
                        }
                
                        // compute fitness for agents
                        
                        fitness[i] = self.problem.objectivefunction(&c[i].genes); //fobj(&c[i]);
                
                        // check fitness with best 
                        if fitness[i] < ceq1_fit {
                            ceq1_index = i;
                            ceq1_fit= fitness[i];
                            copy_vector(&c[i].genes, &mut ceq1, dim);
                        }
                        else if (fitness[i] < ceq2_fit) & (fitness[i] > ceq1_fit) {
                            //ceq2_index = i;
                            ceq2_fit= fitness[i];
                            copy_vector(&c[i].genes, &mut ceq2, dim);            
                        }
                        else if (fitness[i] < ceq3_fit) & (fitness[i] > ceq2_fit) & (fitness[i] > ceq1_fit) {
                            //ceq3_index = i;
                            ceq3_fit= fitness[i];
                            copy_vector(&c[i].genes, &mut ceq3, dim);
                        }
                        else if (fitness[i] < ceq4_fit) & (fitness[i] > ceq3_fit) & (fitness[i] > ceq2_fit) & (fitness[i] > ceq1_fit) {
                            //ceq4_index = i;
                            ceq4_fit= fitness[i];
                            copy_vector(&c[i].genes, &mut ceq4, dim);
                        }
                    }
    
                    //-- Memory saving---
                
                    if iter == 0 {
                        copy_vector(&fitness, &mut fit_old, dim);
                        copy_matrix(&c, &mut c_old);
                    }
                
                    for i in 0..particles_no {
                        if fit_old[i] < fitness[i] {
                            fitness[i]=fit_old[i];
                            copy_vector2genome(&c_old[i], &mut c[i]);
                        }
                    }
                
                    copy_matrix(&c, &mut c_old);
                    copy_vector(&fitness, &mut fit_old, dim);
                
                    // compute averaged candidate Ceq_ave 
                    for i in 0..dim {
                        ceq_ave[i] = (ceq1[i] + ceq2[i] + ceq3[i] + ceq4[i])/4.0;    
                    }
                
                    //Equilibrium pool
                    for i in 0..dim {
                        c_pool[0][i] = ceq1[i];
                        c_pool[1][i] = ceq2[i];
                        c_pool[2][i] = ceq3[i];
                        c_pool[3][i] = ceq4[i];
                        c_pool[4][i] = ceq_ave[i];
                    }
                
                    // comput t using Eq 09
                    let tmpt = (iter / max_iter) as f64;
                    let t : f64 = (1.0 - tmpt).powf(a2*tmpt);
    
                    // let chronos = Instant::now();
                    
                    for i in 0..particles_no {
                
                        randomize(&mut lambda);        //  lambda=rand(1,dim);  lambda in Eq(11)
                        randomize(&mut r);             //  r=rand(1,dim);  r in Eq(11  
                                
                        //-------------------------------------------------------
                        // Ceq=C_pool(randi(size(C_pool,1)),:); 
                        // random selection of one candidate from the pool
                        _index = interval.sample(&mut rng);
                        copy_vector(&c_pool[_index], &mut ceq, dim);
                        //--------------------------------------------------------
                        // compute F using Eq(11) 
                        for j in 0..dim {
                        f[j]=a1*f64::signum(r[j]-0.5)*(f64::exp(-1.0*lambda[j]*t)-1.0); 
                    }
                
                    // r1 and r2 to use them in Eq(15)
                        randomize(&mut r1);
                        randomize(&mut r2);
                
                    for j in 0..dim {
                        // Eq. 15
                        if r2[j]>gp { _gcp =0.5*r1[j]; }
                        else {_gcp =0.0f64;}
                    
                        // Eq. 14
                        _g0 = _gcp*(ceq[j]-lambda[j]*c[i].genes[j]);
                        
                        // Eq 13
                        _g =_g0*f[j];
                        
                        // Eq. 16
                        c[i].genes[j] = ceq[j]+(c[i].genes[j]-ceq[j])*f[j] +  (_g/(lambda[j]*v))*(1.0-f[j]); 
                        }    
                    }
    
                    // let duration = chronos.elapsed();
                    // println!("seq--> End computation in : {:?}", duration);
                
                    convergence_curve[iter] = ceq1_fit;
                    iter+=1;  
    
                    #[cfg(feature = "report")]
                    println!("Iter : {}, Best-fit : {}, Best-solution : {:?}", iter, ceq1_fit, ceq1);
                }
                
                //return results
                let duration = chronos.elapsed();
                let result = OptimizationResult{
                    best_genome : Some(Genome::from(ceq1_index, &ceq1,ceq1_fit)),
                    best_fitness : Some(ceq1_fit),
                    convergence_trend : Some(convergence_curve),
                    computation_time : Some(duration),
                    err_report : None,
                };
                return result;   
            }

        }
    }



}


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
    ub : &'a [f64], a1 :f64, a2 :f64, gp : f64)-> Result<EOparams<'a>, String> {
                      
        let params = EOparams{
            population_size : p_size,
            dimensions : dim,
            max_iterations : max_iter,
            lower_bounds : lb,
            upper_bounds : ub,
            a1,
            a2,
            gp,
        };

       match params.check() {
           Err(error)=> Err(error),
           Ok(())=> Ok(params),
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

impl<'a> Default for EOparams<'a>{

    ///
    /// Return default values of parameters, as following :
    /// 
    /// ~~~
    /// 
    ///  use sefar::sequential_algos::eo::*;
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
    fn default()->Self{
        EOparams{
            population_size : 10,
            dimensions : 3,
            max_iterations : 100,
            lower_bounds : &[-100.0f64, -100.0, -100.0],
            upper_bounds : &[100.0f64, 100.0, 100.0],
            a1 : 2.0f64,
            a2 : 1.0f64,
            gp : 0.5f64,
        }
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

    #[test]
    fn eoparams_unwrap_or_default_test_1(){        
        let _ub =  vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let _lb =  vec![-1.0f64, -2.0, -3.0, -4.0, -5.0];
        
        let p = EOparams::new(10,10,100, _lb.as_slice(), _ub.as_slice(),0.5, 0.5,0.5)
                                       .unwrap_or_default();     
        assert_eq!(p.a1, 2.0f64);
        assert_eq!(p.a2, 1.0f64);
        assert_eq!(p.gp, 0.50f64);
    }

    #[test]
    fn eoparams_unwrap_or_default_test_2(){        
        let _ub =  vec![1.0f64, 2.0, 3.0, 4.0, 5.0];
        let _lb =  vec![-1.0f64, -2.0, -3.0, -4.0, -5.0];
        
        let p = EOparams::new(10,5,100, _lb.as_slice(), _ub.as_slice(),0.5, 0.5,0.5)
                                       .unwrap_or_default();     
        assert_eq!(p.a1, 0.50f64);
        assert_eq!(p.a2, 0.50f64);
        assert_eq!(p.gp, 0.50f64);
    }
}








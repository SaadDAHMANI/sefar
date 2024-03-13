/// -------------------------------------------------------------------------------------
/// Plant competition optimization (PCO)
/// Implemented in Rust programming language by Saad Dahmani (s.dahmani@univ-bouira.dz)
/// Original Matlab code :
/// https://github.com/iman-aliabdi/PCO-Plant-Competition-Optimization
/// Plant competition optimization: A novel metaheuristic algorithm
/// Piper's link : https://onlinelibrary.wiley.com/doi/abs/10.1111/exsy.12956
/// -------------------------------------------------------------------------------------
/// 

extern crate rand;
use rand::distributions::{Distribution, Uniform};

use crate::core::eoa::EOA;
use crate::core::genome::Genome;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
use crate::core::optimization_result::OptimizationResult;
use crate::common::*;

///
/// Plant competition optimization (PCO) (Sequential)
/// Reference:
/// Plant competition optimization: A novel metaheuristic algorithm
/// Piper's link : https://onlinelibrary.wiley.com/doi/abs/10.1111/exsy.12956
/// Original Matlab code:
/// https://github.com/iman-aliabdi/PCO-Plant-Competition-Optimization
/// 
#[derive(Debug)]
pub struct PCO<'a, T : Problem> {
     pub problem : &'a mut T,
     pub params : &'a PCOparams<'a>,
     pub optimization_result : OptimizationResult,
}

impl<'a, T : Problem> PCO<'a, T>{

    pub fn new(settings :&'a PCOparams, problem : &'a mut T )->Self{       
        let result = OptimizationResult{
            best_genome : None,
            best_fitness :None,
            convergence_trend : None,
            computation_time : None,
            err_report : None, 
        };
       
        PCO { 
             problem,
             params: settings,
             optimization_result: result,            
        }
    }
}


impl<'a, T: Problem> EOA for PCO<'a, T> {
   
    fn run(&mut self)-> OptimizationResult{

        let n : usize = self.params.get_population_size();
        let dim : usize = self.params.get_dimensions();
        let ub  = self.params.get_upper_bounds();
        let lb = self.params.get_lower_bounds();

        let noi : usize = self.params.get_max_iterations();
        let max_plant_number : usize = self.params.max_plant_number;
        let vmax : f64 = self.params.vmax as f64;
        let alpha : f64 = self.params.alpha;
        

        

        //----------------------------------------
        let mut rmax_vec : Vec<f64> = vec![0.0f64; dim];
        let mut rmax : f64 = 0.0;
        let mut r : Vec<f64> = vec![0.0f64; n];

        let mut v : Vec<f64> = vec![0.0f64; n];
        let mut best : Vec<f64> = Vec::new();
                
        let mut f : Vec<f64> = vec![0.0f64; n];
        let mut fn_vec : Vec<f64> = vec![0.0f64; n];
        let mut fitness : Vec<f64> = vec![0.0f64; n];
        let mut fc : Vec<f64> = vec![0.0f64; n];
        
        let migrant_seeds_no : usize = 0;
        let migrant_plant : Vec<Genome> = Vec::new();



        let maxteta : f64 = f64::exp(-1.0);
        let teta : f64 = maxteta.clone();
        let mut plant_number= n.clone();
        let mut iteration : usize = 1;

        let max_plant : usize = n.clone();


        //-----------------------------------------
        
        /*  let mut i : usize = 0;
        for (a,b) in ub.iter().zip(lb.iter()){
            rmax[i] = a-b;
            i+=1;
        } */

        rmax = ub[0]-lb[0];

        println!("rmax : {:?}", rmax);

        //-------------------------------------------


        let mut plants = self.initialize(self.params);
        let x= plants.clone();

        randomize(&mut v);

        while plant_number <= max_plant_number && iteration <= noi {

            //for i=1:plantNumber
                //f(i)=fobj(plants(i,:));
            //end

            // Evaluation of candidate solutions:
            for i in 0..plants.len() {
                f[i] = self.problem.objectivefunction(&plants[i].genes);
            }

            // Calculate Fitness Coefficient=fc
            let minf = f.iter().fold(f64::MAX, |minf, y| minf.min(*y));
            best.push(minf);

            let normf = f.iter().map(|&x| x * x).sum::<f64>().sqrt();

            for i in 0..n {
                fn_vec[i] = f[i]/normf;                
            }

            for i in 0..n {
                fitness[i] = 1.0/(1.0 + fn_vec[i]);
            }

            let mx : f64 = fitness.iter().fold(f64::MIN, |mx, y| mx.max(*y));
            let mn : f64 = fitness.iter().fold(f64::MAX, |mn, y| mn.min(*y));

            //println!("mx : {}; mn : {}", mx,mn);
            if mx == mn {
                for i in 0..n {
                    fc[i] = fitness[i]/mx;
                }    
            }
            else {
                // fc=(fitness-mn)./(mx-mn);
                let dif_mx_mn = mx-mn;
                for i in 0..n {
                    fc[i] = (fitness[i]-mn)/ dif_mx_mn;
                }
            }
            
            fitness.sort_by(|a,b| b.partial_cmp(a).unwrap());

            #[cfg(feature="report")] println!("fitness : {:?}", fitness);

            let mut survive : Vec<bool> = Vec::new();
            for &value in &fc {
                survive.push(value >= fc[max_plant-1]);
            }

            #[cfg(feature="report")] println!("survive : {:?}", survive );
          

            let mut  new_plant : Vec<Genome>= Vec::new();
            
            for i in 0..plants.len() {
                if survive[i] == true {
                    new_plant.push(plants[i].clone());
                }
            }

            //plants=newPlant;
            plants = new_plant;

            //sz=size(newPlant);
            //let sz : usize = plants.len();
            //let mut x1 : Vec<f64> = Vec::new();
            //let mut y : Vec<f64> = Vec::new();

            //for j in 0..plants.len() {
                //x1.push(plants[j].genes[0]);
                //y.push(plants[j].genes[1]);
            //}

            plant_number = plants.len(); //x1.len();

            #[cfg(feature ="report")] println!("Iter {};  plant_number :{}", iteration, plant_number);

            // st=zeros(plantNumber,1);   
            let mut st : Vec<f64> = vec![0.0f64; plant_number];


            for i in 0..plant_number {
                //Compute Neighborhood Radius
                //r(i)=teta*rmax*exp(1-(5*v(i))/vmax);

                r[i] = teta * rmax*f64::exp(1.0-(alpha*v[i])/vmax);

                
                // non : number of neighbours
                let mut non : usize =0;

                for j in 0..plant_number {
                    if euclidian_dist(&plants[i], &plants[j]) <= r[i] { // are neighbours in this case:
                        // st(i)=st(i)+v(j);
                        //non=non+1;
                        st[i] += v[j];
                        non +=1;                  
                    }
                }

                

                  


              
                println!("r : {:?}", r);                
            }




              //plant_number += 1;
            iteration +=1; 
        }



        let result = OptimizationResult{
            best_genome : None, //Some(Genome::new(0, self.params.dimensions)),
            best_fitness : None, //Some(-1111.1111),
            convergence_trend : None, //Some(convergence_curve),
            computation_time : None, //Some(duration),
            err_report : None,
        };
        return result;   

    }    

}






#[derive(Debug, Clone)]
pub struct PCOparams<'a> {
    /// Number of initial plants.
    pub population_size : usize,

    /// Size of the problem (number of decision variables).
    pub dimensions : usize,

    /// Maximum number of iterations.
    pub max_iterations: usize,

    pub lower_bounds : &'a [f64],
    pub upper_bounds : &'a [f64],
    
    /// vmax : Maximum size of plants.
    pub vmax : usize,

    /// Maximum of Plants Number.
    pub max_plant_number : usize,

    /// Growth rate.
    /// It shoud be greather than the maximum value, Theta_max = e^-1=0.36788.
    pub theta : f64,


    /// Parameter (can be considered as constant too).
    /// Defautl value : alpha = 5.0. 
    pub alpha : f64,
    
    /// Growth parameter. 
    pub k : f64,

    ///Seed migration rate.
    pub miu : f64,    
}

#[allow(dead_code)]
impl<'a> PCOparams<'a>{
    pub fn new(p_size: usize, dim : usize, max_iter : usize, lb : &'a [f64], 
    ub : &'a [f64], vmax : usize, max_plant_number : usize, theta : f64, alpha : f64, k : f64, miu : f64)-> Result<PCOparams<'a>, String> {
                      
        let params = PCOparams{
            population_size : p_size,
            dimensions : dim,
            max_iterations : max_iter,
            lower_bounds : lb,
            upper_bounds : ub,
            vmax,
            max_plant_number,
            theta : f64::min(theta, f64::exp(-1.0)),
            alpha,
            k,
            miu,
        };

       match params.check() {
           Err(error)=> Err(error),
           Ok(())=> Ok(params),
       }
    }
}

impl<'a> Parameters for PCOparams<'a> {

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
 

impl<'a> Default for PCOparams<'a>{

    ///
    /// Return default values of parameters, as following :
    /// 
    /// ~~~
    /// 
    ///  use sefar::sequential_algos::pco::*;
    /// 
    ///  PCOparams{
    ///     population_size : 10,
    ///     dimensions : 3,
    ///     max_iterations : 100,
    ///     lower_bounds : &[-100.0f64, -100.0, -100.0],
    ///     upper_bounds : &[100.0f64, 100.0, 100.0],
    ///     vmax : 20,
    ///     max_plant_number : 500,
    ///     theta : 0.005,
    ///     alpha : 5.0,
    ///     k : 0.1,
    ///     miu : 0.05,
    /// };
    /// ~~~
    /// 
    fn default()->Self{
        PCOparams{
            population_size : 10,
            dimensions : 3,
            max_iterations : 100,
            lower_bounds : &[-100.0f64, -100.0, -100.0],
            upper_bounds : &[100.0f64, 100.0, 100.0],
            vmax : 20,
            max_plant_number : 500,
            theta : 0.005,
            alpha : 5.0,
            k : 0.1,
            miu : 0.05,
        }
    }
}


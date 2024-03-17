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
use rand::Rng;

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
#[allow(dead_code)]
pub struct PCO<'a, T : Problem> {
     pub problem : &'a mut T,
     pub params : &'a PCOparams<'a>,
     pub optimization_result : OptimizationResult,
}

#[allow(dead_code)]
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

#[allow(dead_code)]
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
        let k :f64  = self.params.k;
        let miu : f64 = self.params.miu;

        //----------------------------------------
        //let rmax_vec : Vec<f64> = vec![0.0f64; dim];
        let mut rmax : f64 = 0.0;
        let mut r : Vec<f64> = Vec::new();

        let mut v : Vec<f64> = vec![0.0f64; n];
        let mut dv : Vec<f64> = Vec::new();// vec![0.0f64; n];
       
        let mut best : Vec<f64> = Vec::new();

        let mut nos : Vec<usize> = Vec::new();
                
        let mut f : Vec<f64> = Vec::new(); //let mut f : Vec<f64> = vec![0.0f64; n];
        let mut fn_vec : Vec<f64> = Vec::new();

        let mut fitness : Vec<f64> = Vec::new(); //vec![0.0f64; n];
        let mut fc : Vec<f64> = Vec::new(); //vec![0.0f64; n];
        
        let mut migrant_seeds_no : usize = 0;
        let mut migrant_plant :Vec<usize> = Vec::new();

        //--------------------------------------------
        let between = Uniform::from(0.0..=1.0);
        let mut rng = rand::thread_rng();
        //--------------------------------------------

        let maxteta : f64 = f64::exp(-1.0);
        let teta : f64 = maxteta.clone();
        let mut plant_number= n.clone();
        let mut iteration : usize = 1;

        let mut max_plant : usize = n.clone();
        let mut plant_number_old : usize =0;


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

            #[cfg(feature="report")]println!("---- Iter : [{}]", iteration);
            //for i=1:plantNumber
                //f(i)=fobj(plants(i,:));
            //end

            // Evaluation of candidate solutions:
            f.clear();

            for i in 0..plants.len() {
                f.push(self.problem.objectivefunction(&plants[i].genes));
            }

            // Calculate Fitness Coefficient=fc

           let min_f = match f.iter().min_by(|a, b| a.total_cmp(b)){
                Some(min_value) => *min_value,
                None => f64::MAX,
           };

            #[cfg(feature="report")] println!("fitness before sorting, f : {:?}", f);
            #[cfg(feature= "report")] println!("minf = min(f) = {}", min_f);

            best.push(min_f);

            let normf = f.iter().map(|&x| x * x).sum::<f64>().sqrt();

           fn_vec.clear();
            for i in 0..f.len() {
                fn_vec.push(f[i]/normf);                
            }

            fitness.clear();

            for i in 0..fn_vec.len() {
                fitness.push(1.0/(1.0 + fn_vec[i]));
            }

            let mx : f64 = fitness.iter().fold(f64::MIN, |mx, y| mx.max(*y));
            let mn : f64 = fitness.iter().fold(f64::MAX, |mn, y| mn.min(*y));

            //println!("mx : {}; mn : {}", mx,mn);
            
            fc.clear();
            
            if mx == mn {
                for i in 0..fitness.len() {
                    fc.push(fitness[i]/mx);
                }    
            }
            else {
                // fc=(fitness-mn)./(mx-mn);
                let dif_mx_mn = mx-mn;
                for i in 0..fitness.len() {
                    fc.push((fitness[i]-mn)/ dif_mx_mn);
                }
            }
            
            fitness.sort_by(|a,b| b.partial_cmp(a).unwrap());

            #[cfg(feature="report")] println!("fitness after sorting : {:?}", fitness);

            max_plant = plant_number;

            let mut survive : Vec<bool> = Vec::new();
            for &value in &fc {
                survive.push(value >= fc[max_plant-1]);
            }

            #[cfg(feature="report")] println!("survive : {:?}", survive );
          

            let mut  new_plant : Vec<Genome>= Vec::new();
            
            for i in 0..survive.len() {
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

            r.clear();
            dv.clear();

            for i in 0..plant_number {
                //Compute Neighborhood Radius
                //r(i)=teta*rmax*exp(1-(5*v(i))/vmax);

                r.push(teta * rmax*f64::exp(1.0-(alpha*v[i])/vmax));
                
                // non : number of neighbours
                let mut non : f64 =0.0;

                for j in 0..plant_number {
                    if euclidian_dist(&plants[i], &plants[j]) <= r[i] { // are neighbours in this case:
                        // st(i)=st(i)+v(j);
                        //non=non+1;
                        st[i] += v[j];
                        non +=1.0;                  
                    }
                }

                // dv(i)=fc(i)*k*(log(non*vmax)-log(st(i)));
                let tmpdvi = fc[i]*k*(f64::ln(non*vmax)- f64::ln(st[i]));
                dv.push(tmpdvi);

                // if v(i)+dv(i)<vmax
                //      v(i)=v(i)+dv(i);
                // else
                //      v(i)=vmax;
                // end

                if (v[i]+dv[i]) < vmax {
                    v[i] = v[i] + dv[i];
                } 
                else {
                    v[i] = vmax;
                }
            }

             // ----- SEED PRODUCTION ----------------
             // sumNos=0;
             let mut sum_nos : usize =0;
             nos.clear();

             for i in 0..plant_number {
                // NOS(i)=floor(v(i)+1);
                nos.push((v[i]+1.0).floor() as usize);
                
                //sumNos=sumNos+NOS(i);
                sum_nos += nos[i];

                for j in 0..nos[i]{
                    //RND=randi(dim);
                    let rnd = rng.gen_range(0..dim);
                    
                    //Temp=(plants(i,RND)-r(i))+2*r(i)*rand;
                    let rand01 = between.sample(&mut rng);
                    let  temp = (plants[i].genes[rnd]-r[i])+(2.0*r[i]* rand01);
                    
                    // seed=plants(i,:);
                    // seed(RND)=Temp;
                    // plants=[plants;seed];
                    // v=[v;rand];
                    let mut seed = plants[i].clone();
                    seed.genes[rnd]=temp;
                    plants.push(seed);

                    let rand010 = between.sample(&mut rng);
                    v.push(rand010);

                }
             }

             // -- SEED MIGRATION -------------------
             
             // migrantSeedsNoOld = migrantSeedsNo;
             let migrant_seeds_no_old = migrant_seeds_no.clone();
             
             //migrantSeedsNo=floor(miu*sumNos);
             migrant_seeds_no = (miu*sum_nos as f64).floor() as usize;

             //migrantPlantOld=migrantPlant;
             //let migrant_plant_old = migrant_plant.clone();
             if (plant_number+1) <= (plant_number + sum_nos) {

            
                //migrantPlant=randi([plantNumber+1,plantNumber+sumNos],1,migrantSeedsNo);
                migrant_plant = randi(plant_number+1,plant_number + sum_nos, migrant_seeds_no);

                /* for i=1:migrantSeedsNo
                    temp=A+(B-A).*rand(1,dim);
                    plants(migrantPlant(i),:)=temp;
                end */
                
                for i in 0..migrant_seeds_no {
                    for j in 0..dim {
                        plants[migrant_plant[i]].genes[j] = lb[j] + (ub[j]-lb[j])*between.sample(&mut rng); 
                    }
                } 

            }

            //plantNumberOld=plantNumber;
            plant_number_old = plant_number;
            iteration +=1; 
        }

        let result = OptimizationResult{
            best_genome : None, //Some(Genome::new(0, self.params.dimensions)),
            best_fitness : None, //Some(-1111.1111),
            convergence_trend : Some(best), //Some(convergence_curve),
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


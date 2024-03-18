
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
/// QAGO
/// Reference:
/// 
/// 
/// 
/// 
/// 

#[derive(Debug)]
pub struct QAGO<'a, T : Problem> {
     pub problem : &'a mut T,
     pub params : &'a QAGOparams<'a>,
     pub optimization_result : OptimizationResult,
}

impl<'a, T : Problem> QAGO<'a, T>{

    pub fn new(settings :&'a QAGOparams, problem : &'a mut T )->Self{       
        let result = OptimizationResult{
            best_genome : None,
            best_fitness :None,
            convergence_trend : None,
            computation_time : None,
            err_report : None, 
        };
       
        QAGO { 
             problem,
             params: settings,
             optimization_result: result,            
        }
    }

    fn get_empty_solutions(&self, n : usize)->Vec<Genome>{
        let mut result : Vec<Genome> = Vec::with_capacity(n);
        for i in 0..n {
            result.push(Genome::new(i, self.params.get_dimensions()));
        }
        result
    }

    fn select_id(&self, n: usize)->(Vec<usize>, Vec<usize>, Vec<usize>, Vec<usize>){

        let mut l1 : Vec<usize> = vec![0; n];
        let mut l2 : Vec<usize> = vec![0; n];
        let mut l3 : Vec<usize> = vec![0; n];
        let mut l4 : Vec<usize> = vec![0; n];

        let mut r =[0, 0, 0, 0];       

        for i in 0..n{
            if n >= 1 {
                let vecc: Vec<usize> = (0..n).filter(|&x| x != i).collect();
                for kkk in 0..4{
                    let nn = n-kkk;
                    let interval01 = Uniform::from(0..nn);
                    let t = interval01.sample(&mut rand::thread_rng());
                    r[kkk] = vecc[t];
                }

                l1[i] = r[0];
                l2[i] = r[1];
                l3[i] = r[2];
                l4[i] = r[3];
            }
        }
        (l1, l2, l3, l4)
    }   
}

impl <'a, T : Problem> EOA for QAGO<'a, T>{
    fn run(&mut self)-> OptimizationResult {

        let d : usize =  self.params.get_dimensions();
        let n : usize = self.params.get_population_size();
        let max_iter = self.params.get_max_iterations();
        let ub = self.params.get_upper_bounds();
        let lb = self.params.get_lower_bounds();

        let mut gbestx : Genome = Genome::new(n+1, d);      
        let mut gbestfitness : f64 = f64::MAX; // Change this to f64::MIN in the case of maximization.
        let mut gbesthistory = vec![f64::NAN; max_iter];
        let mut fitness = vec![0.0; n];
        let mut iter : usize =0;

        let mut p2 : Vec<f64> = vec![0.0; n];
        let mut p3 : Vec<Vec<f64>> = vec![vec![0.0; d]; n];

        let n_f64 = n as f64;
        let between01 = Uniform::from(0.0..=1.0);

        let mut best_x : Genome = Genome::new(n+2, d);
        let mut gap : Vec<Vec<f64>> = vec![vec![0.0; d]; 5];

        let mut worst_x : Vec<Genome> = self.get_empty_solutions(n);
        let mut better_x : Vec<Genome> = self.get_empty_solutions(n);
        let mut normal_x : Vec<Genome> = self.get_empty_solutions(n);

        //let mut better_x : Vec<Genome> = Vec::new();

        // Step 1 : Initialization
        let mut x = self.initialize(self.params);

        //Evaluation of candidate solution using the objective function
        for i in 0..n{
            fitness[i] = self.problem.objectivefunction(&x[i].genes);
            
            if gbestfitness >= fitness[i] {
                gbestfitness = fitness[i];
                copy_vector(&x[i].genes, &mut gbestx.genes, d);                
            }            
        }
        // Save the best fitness history:
        gbesthistory[iter]= gbestfitness;


        // Loop iterations
        while iter < max_iter {

           //Sorte and sorting indexes:
           let mut ind : Vec<usize> = (0..fitness.len()).collect();
           ind.sort_by(|&a, &b| fitness[a].partial_cmp(&fitness[b]).unwrap());
           //------------------------------
           // Parameter adaptation based on distribution
           // P1=ceil(unifrnd(0.05,0.2)*N);
           let p1 = f64::ceil(uniform_rand1(0.05, 0.2)*n as f64);
           let p1_usize = p1.round() as usize;

           #[cfg(feature="report")] println!("p1 = {}", p1);
           // P2=normrnd(0.001*ones(1,N),0.001*ones(1,N));
            for j in 0..n{
                p2[j] = normal_rand1(0.001*n_f64, 0.001*n_f64);
            }

            #[cfg(feature="report")] println!("QAGO : P2 = {:?}", p2);

            //P3=normrnd(0.3*rand(N,D),0.01);
            for i in 0..n{
                for j in 0..d {
                    p3[i][j] = normal_rand1(0.3*between01.sample(&mut rand::thread_rng()), 0.01)
                }
            }
           
            #[cfg(feature="report")] println!("QAGO : P3 = {:?}", p3);

            //1. Improved learning phase 
            //1.1 Sampling individuals 

            //Best_X=x(ind(1),:);
             copy_vector(&x[ind[0]].genes, &mut best_x.genes, d); 
            //-------------------------------------------------------------------------------------
             //worse_index=ind(randi([N-P1+1,N],N,1));
            let worse_index = rand_vec(n-p1_usize+1, n, n);
                        
            //Worst_X=x(worse_index,:);            
             for k in 0..n {
                 copy_vector(&x[worse_index[k]].genes, &mut worst_x[k].genes, d);//  worst_x[*k] = (x[worse_index[*k]].clone());
             }
           
             println!("_______________________________________________________________________");

             #[cfg(feature = "report")] println!("Worst_index : {:?}; Worst_X : {:?}", worse_index, worst_x);

             //-------------------------------------------------------------------------------------
             //better_index=ind(randi([2,P1],N,1));

             let better_index = rand_vec(2, p1_usize, n);
             
             //Better_X=x(better_index,:);
             for k in 0..n {
                 copy_vector(&x[better_index[k]].genes, &mut better_x[k].genes, d);
             }

             #[cfg(feature="report")] println!("better_x : {:?}", better_x);            
            //-------------------------------------------------------------------------------------
            //normal_index=ind(randi([P1+1,N-P1],N,1));
            let normal_index = rand_vec(p1_usize+1, n-p1_usize, n);
            
            //Normal_X=x(normal_index,:);
            for k in 0..n {
                 copy_vector(&x[normal_index[k]].genes, &mut normal_x[k].genes, d); 
            };
            #[cfg(feature="report")] println!("normal_x : {:?}", normal_x);
            //------------------------------------------------------------------------------------- 

            //[L1,L2,L3,L4]=selectID(N);
            let (l1, l2, l3, l4) = self.select_id(n);

            println!("_______________________________________________________");

            #[cfg(feature="report")] println!("l1 : {:?}", l1);

            for i in 0..n {
                for j in 0..d {
                    //  Gap(1,:)=(Best_X - Better_X(i,:));
                    gap[0][j] = best_x.genes[j]-better_x[i].genes[j];

                    //Gap(2,:)=(Better_X(i,:)-Normal_X(i,:));
                    gap[1][j] = better_x[i].genes[j] - normal_x[i].genes[j];

                    // Gap(3,:)=(Normal_X(i,:)-Worst_X(i,:));
                    gap[2][j] = normal_x[i].genes[j] - worst_x[i].genes[j];

                    //Gap(4,:)=(x(L1(i),:)-x(L2(i),:));
                    gap[3][j] = x[l1[i]].genes[j] - x[l2[i]].genes[j];

                    //Gap(5,:)=(x(L3(i),:)-x(L4(i),:));
                    gap[4][j] = x[l3[i]].genes[j] - x[l4[i]].genes[j];
                }

                // Parameter self-adaptation based on one-dimensional mapping of vectors

               

            }


            





            //iteration incrementation
            iter+=1;
        }





        let result = OptimizationResult{
            best_genome : None,
            best_fitness : None, 
            convergence_trend : None,
            computation_time : None,
            err_report : None,
        };
        return result;  


    }


}






#[derive(Debug, Clone)]
pub struct QAGOparams<'a> {
    pub population_size : usize,
    pub dimensions : usize,
    pub max_iterations: usize,
    pub lower_bounds : &'a [f64],
    pub upper_bounds : &'a [f64],
}

impl<'a> Parameters for QAGOparams<'a> {}

impl<'a> Default for QAGOparams<'a>{

    ///
    /// Return default values of parameters, as following :
    /// 
    /// ~~~
    /// 
    ///  use sefar::sequential_algos::qago::*;
    /// 
    ///  QAGOparams{
    ///     population_size : 10,
    ///     dimensions : 3,
    ///     max_iterations : 100,
    ///     lower_bounds : &[100.0f64, 100.0, 100.0],
    ///     upper_bounds : &[-100.0f64, -100.0, -100.0],
    /// };
    /// ~~~
    /// 
    fn default()->Self{
        QAGOparams {
            population_size : 4,
            dimensions : 3,
            max_iterations : 1,
            lower_bounds : &[-100.0f64, -100.0, -100.0],
            upper_bounds : &[100.0f64, 100.0, 100.0],
        }
    }
}

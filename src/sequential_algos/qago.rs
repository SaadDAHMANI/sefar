
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
        //let mut worst_x : Vec<Genome> = Genome::new(n+2, d);
        //let mut better_x : Genome = Genome::new(n+2, d);

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

             //worse_index=ind(randi([N-P1+1,N],N,1));
             let tmp_worse_index = rand_vec(n-p1_usize+1, n, n);

             let mut worse_index : Vec<usize> = Vec::new();
             for k in tmp_worse_index {
                worse_index.push(ind[k]);
             }

             //Worst_X=x(worse_index,:);

             let mut worst_x : Vec<Genome> = Vec::new();
             for k in worse_index.iter() {
                worst_x.push(x[worse_index[*k]].clone());
             }
           
             println!("_______________________________________________________________________");

             #[cfg(feature = "report")] println!("Worst_index : {}; Worst_X : {:?}", worse_index, worst_x);

             //better_index=ind(randi([2,P1],N,1));

             let tmp_better_index = rand_vec(2, p1_usize, n, );
             let mut better_index : Vec<usize> = Vec::new();
             for k in tmp_better_index {
                better_index.push(ind[k]);
             }

              //Better_X=x(better_index,:);
             let mut better_x : Vec<Genome> = Vec::new();
             for k in better_index.iter() {
                better_x.push(x[better_index[*k]].clone());
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

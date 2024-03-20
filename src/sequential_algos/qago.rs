
extern crate rand;
use rand::distributions::{Distribution, Uniform};
use rand::rngs::ThreadRng;
use rand::Rng;
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

    fn select_id(&self, index_differ : usize, n : usize, rng : &mut ThreadRng)->(usize, usize, usize, usize){

        let mut l1 : usize = index_differ; 
        let mut l2 : usize = index_differ;
        let mut l3 : usize = index_differ;
        let mut l4 : usize = index_differ;

        //let mut rng = rand::thread_rng();
        let mut interval = Uniform::from(0..n);        
        
        while l1 == index_differ {
            l1 = interval.sample(rng);
        }

        while l2 == index_differ {
            l2 = interval.sample(rng);
        }

        while l3 == index_differ {
            l3 = interval.sample(rng);
        }

        while l4 == index_differ {
            l4 = interval.sample(rng);
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
        let mut rng = rand::thread_rng();

        let mut best_x : Genome = Genome::new(n+2, d);
        let mut gap : Vec<Vec<f64>> = vec![vec![0.0; d]; 5];
        let mut dgap : [f64; 5] = [0.0; 5];
        let mut fgap : [f64; 5] = [0.0; 5];

        let mut lf : [f64; 5] = [0.0; 5];
        let mut sf: [f64; 5]  = [0.0; 5];
        let mut ls : [f64; 5] = [0.0; 5];
       
        let mut learn_operator : Vec<f64> = vec![0.0; d];
       
        let mut worst_x : Vec<Genome> = self.get_empty_solutions(n);
        let mut better_x : Vec<Genome> = self.get_empty_solutions(n);
        let mut normal_x : Vec<Genome> = self.get_empty_solutions(n);
        let mut newx : Vec<Genome> = self.get_empty_solutions(n);

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
            //let worse_index = rand_vec(n-p1_usize+1, n-1, n);
            let worse_index = rand_vec(n-p1_usize, n-1, n);
            
            println!("worse_index : {:?}", worse_index);
                        
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
         
             println!("_______________________________________________________");

            for i in 0..n {
                  //[L1,L2,L3,L4]=selectID(N);

                let (l1, l2, l3, l4) = self.select_id(i, n , &mut rng);
                
                #[cfg(feature="report")] println!("i : {}, l1: {}, l2: {}, l3: {}, l4: {}", i, l1, l2, l3, l4);

                for j in 0..d {
                    //  Gap(1,:)=(Best_X - Better_X(i,:));
                    gap[0][j] = best_x.genes[j]-better_x[i].genes[j];

                    //Gap(2,:)=(Better_X(i,:)-Normal_X(i,:));
                    gap[1][j] = better_x[i].genes[j] - normal_x[i].genes[j];

                    // Gap(3,:)=(Normal_X(i,:)-Worst_X(i,:));
                    gap[2][j] = normal_x[i].genes[j] - worst_x[i].genes[j];

                    //Gap(4,:)=(x(L1(i),:)-x(L2(i),:));
                    gap[3][j] = x[l1].genes[j] - x[l2].genes[j];

                    //Gap(5,:)=(x(L3(i),:)-x(L4(i),:));
                    gap[4][j] = x[l3].genes[j] - x[l4].genes[j];
                }

                // Parameter self-adaptation based on one-dimensional mapping of vectors

                dgap[0] = best_x.genes.iter().zip(better_x[i].genes.iter()).fold(0.0f64, |sum, (a, b)| sum + (a*b));
                dgap[1] = better_x[i].genes.iter().zip(normal_x[i].genes.iter()).fold(0.0f64, |sum, (a, b)| sum + (a*b));
                dgap[2] = normal_x[i].genes.iter().zip(worst_x[i].genes.iter()).fold(0.0f64, |sum, (a, b)| sum + (a*b));
                dgap[3] = x[l1].genes.iter().zip(x[l2].genes.iter()).fold(0.0f64, |sum, (a, b)| sum + (a*b));
                dgap[4] = x[l3].genes.iter().zip(x[l4].genes.iter()).fold(0.0f64, |sum, (a, b)| sum + (a*b));

                let min_distance : f64 = match dgap.iter().min_by(|a, b| a.partial_cmp(b).unwrap()){
                    Some(value) => (*value*2.0).abs(),
                    None => 1.0,
                };

                //DGap=DGap+2*abs(minDistance)
                //let min_distance_2 =  2.0*min_distance.abs();

                for j in 0..5 {
                    dgap[j] +=  min_distance;
                }

                let sum_dgap = dgap.iter().fold(0.0f64, |sum, a| sum+a); 
                for k in 0..5{
                    lf[k]= dgap[k]/sum_dgap + 1.0;
                }

                //Parameter self-adaptation based on fitness difference
                //FGap(1,:)=(abs(fitness(ind(1))-fitness(better_index(i))));
                fgap[0] = (fitness[ind[0]]-fitness[better_index[i]]).abs();

                //FGap(2,:)=(abs(fitness(better_index(i))-fitness(normal_index(i))));
                fgap[1] = (fitness[better_index[i]] - fitness[normal_index[i]]).abs();

                //FGap(3,:)=(abs(normal_index(i)-fitness(worse_index(i)))); // Err in the original code line.
                fgap[2] = (fitness[normal_index[i]] - fitness[worse_index[i]]).abs();

                //FGap(4,:)=(abs(fitness(L1(i))-fitness(L2(i))));
                fgap[3] = (fitness[l1] - fitness[l2]).abs();
                
                //FGap(5,:)=(abs(fitness(L3(i))-fitness(L4(i))));
                fgap[4] = (fitness[l3] - fitness[l4]).abs();

                //SF=FGap./sum(FGap);
                let sum_fgap = fgap.iter().fold(0.0f64, |sum, a| sum + a);
                for k in 0..5 {
                    sf[k] = fgap[k]/sum_fgap+1.0;
                }

                //Parameter self-adaptation based on Jensen-Shannon divergence
                // LS=(LF+SF)/2;
                for k in 0..5 {
                    ls[k] = (sf[k] + lf[k])/2.0;
                }

                //Djs=0.5*sum(LF.*log(LF./LS))+0.5*sum(SF.*log(SF./LS));

                let sum1 : f64 = lf.iter().zip(ls.iter()).fold(0.0f64, |sum, (lfa,lsb)| sum + (lfa*f64::ln(lfa/lsb)));
                let sum2 : f64 = sf.iter().zip(lf.iter()).fold(0.0f64, |sum, (sfa, lfb)| sum + (sfa*f64::ln(sfa/lfb)));

                // djs=sqrt(Djs);
                let djs : f64 = f64::sqrt(0.5*(sum1 + sum2)); 

                //Learning operator refinement
                //newx(i,:)=x(i,:)+sum(Gap.*(djs.*LF+(1-djs).*SF),1);
                for j in 0..5 {
                    let multiplier: f64 = djs*lf[j] + (1.0-djs)*sf[j];
                    for l in 0..d {
                        gap[j][l] *= multiplier; 
                    }
                }

                 for t in 0..d {
                    let mut tmp_sum : f64 =0.0;
                    for j in 0..5 {
                        tmp_sum += gap[j][t];
                    }
                    learn_operator [t] = tmp_sum;
                }   

                #[cfg(feature = "report")] println!("lear_operator : {:?}", learn_operator);
                                      











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

impl<'a> Parameters for QAGOparams<'a> {
    fn get_dimensions(&self)->usize {
         self.dimensions    
    }

    fn get_max_iterations(&self)->usize {
        self.max_iterations
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
            population_size : 10,
            dimensions : 3,
            max_iterations : 1,
            lower_bounds : &[-100.0f64, -100.0, -100.0],
            upper_bounds : &[100.0f64, 100.0, 100.0],
        }
    }
}

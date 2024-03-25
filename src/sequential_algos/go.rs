
use crate::core::eoa::EOA;
use crate::core::genome::Genome;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
use crate::core::optimization_result::OptimizationResult;
use crate::common::*;

///
/// GO : Growth Optimizer  
/// Reference:
/// 
/// 
/// 
/// 
/// 
#[derive(Debug)]
pub struct GO<'a, T : Problem> {
     pub problem : &'a mut T,
     pub params : &'a GOparams<'a>,
     pub optimization_result : OptimizationResult,
}


impl<'a, T : Problem> GO<'a, T> {

    pub fn new(settings :&'a GOparams, problem : &'a mut T )->Self{       
        let result = OptimizationResult{
            best_genome : None,
            best_fitness :None,
            convergence_trend : None,
            computation_time : None,
            err_report : None, 
        };
       
        GO{ 
             problem,
             params: settings,
             optimization_result: result,            
        }
    }
}

impl<'a, T : Problem> EOA for GO<'a, T> {
    fn run(&mut self)-> OptimizationResult {

        let n : usize = self.params.population_size;
        let d : usize = self.params.dimensions;
        let max_iter : usize = self.params.max_iterations;

        let mut iter : usize = 0;        
        //Parameter setting
        const P1 : usize = 5;
        const P2 :f64 = 0.001;
        const P3 :f64 = 0.3;

        let mut fes : usize = 0;
        let mut gbestfitness : f64 = f64::MAX;

        let mut fitness : Vec<f64> = vec![0.0; n];
        let mut gbest_x : Genome = Genome::new(n+1, d);
        let mut gbesthistory : Vec<f64> = vec![0.0; max_iter];

        let mut best_x : Genome = Genome::new(n+2, d); 


        //Initialization
        let mut x = self.initialize(self.params);
        

        //Evaluation of search agents

        for i in 0..n {
            fitness[i] = self.problem.objectivefunction(&x[i].genes);
            fes +=1;

            if gbestfitness > fitness[i] {
                gbestfitness = fitness[i];
                copy_vector(&x[i].genes, &mut gbest_x.genes,d);
            }
        }

        gbesthistory[0] = gbestfitness;

        println!("Best_fitness : {}", gbestfitness);


        while  iter < max_iter {

            //Sorte and sorting indexes:
           let mut ind : Vec<usize> = (0..n).collect();
           ind.sort_by(|&a, &b| fitness[a].total_cmp(&fitness[b]));
           //----------------------------------------------------------------------------------------

           // Save best solution 
           copy_vector(&x[ind[0]].genes, &mut best_x.genes, d);

            // Learning phase

            for i in 0..n {
                let worse_index = rand_vec(n-P1, n-1, n);

                println!("wors_index : {:?}", worse_index);


            }


            


            iter += 1;
        }










        let result : OptimizationResult = OptimizationResult {
            best_genome : None,
            best_fitness : None,
            convergence_trend : None,
            computation_time : None,
            err_report : None,
        };

        result

    }
} 




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
        usize::max(self.population_size, 6)
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


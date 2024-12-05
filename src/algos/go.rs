use crate::common::*;
use crate::core::eoa::{InitializationMode, EOA};
use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
use rand::rngs::ThreadRng;
use rand_distr::{Distribution, Uniform};
#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use std::time::Instant;

///
/// GO : Growth Optimizer & Binary-GO
/// Reference:
/// Growth Optimizer: A powerful metaheuristic algorithm for solving continuous and
/// discrete global optimization problems (2023).
/// Paper link : <https://doi.org/10.1016/j.knosys.2022.110206>
/// Matlab original code : <https://github.com/tsingke/Growth-Optimizer>
///
/// Written in Rust by Saad Dahmani <sd.dahmani2000@gmail.com>
///
#[derive(Debug)]
pub struct GO<'a, T: Problem> {
    /// The problem to optimize. It must define the Problem trait.
    pub problem: &'a mut T,

    /// Define the parameters of GO algorithm.
    pub params: &'a GOparams<'a>,
}

impl<'a, T: Problem> GO<'a, T> {
    ///
    /// Return a new instance of the Growth Optimizer (GO) algorithm.
    /// settings: The optimization parameters,
    /// problem: The problem to optimize.
    ///
    pub fn new(settings: &'a GOparams, problem: &'a mut T) -> Self {
        GO {
            problem,
            params: settings,
        }
    }

    fn get_empty_solutions(&self, n: usize) -> Vec<Genome> {
        let mut result: Vec<Genome> = Vec::with_capacity(n);
        for i in 0..n {
            result.push(Genome::new(i, self.params.get_problem_dimension()));
        }
        result
    }

    ///
    /// Generate 2 random values in [0, maxi[ differ from index_differ.
    ///
    fn select_id(&self, index_differ: usize, maxi: usize, rng: &mut ThreadRng) -> (usize, usize) {
        let mut l1: usize = index_differ;
        let mut l2: usize = index_differ;

        //let mut rng = rand::thread_rng();
        let interval = Uniform::from(0..maxi);

        while l1 == index_differ {
            l1 = interval.sample(rng);
        }

        while l2 == index_differ {
            l2 = interval.sample(rng);
        }

        (l1, l2)
    }

    fn norm(&self, a: &[f64]) -> f64 {
        a.iter().fold(0.0, |sum, x| sum + (x * x)).sqrt()
    }
}

impl<'a, T: Problem> EOA for GO<'a, T> {
    ///
    /// Call this function to execute GO algorithm.
    ///
    fn run(&mut self) -> OptimizationResult {
        let chronos = Instant::now();

        let n: usize = self.params.population_size;
        let d: usize = self.params.problem_dimension;
        let max_iter: usize = self.params.max_iterations;

        let ub = self.params.upper_bounds;
        let lb = self.params.lower_bounds;

        let mut iter: usize = 0;
        //Parameter setting
        const P1: usize = 5;
        const P2: f64 = 0.001;
        const P3: f64 = 0.3;

        let mut fes: f64 = 0.0;
        let max_fes: f64 = n as f64 + (max_iter * 2 * n) as f64;

        let mut gbestfitness: f64 = f64::MAX;

        let mut fitness: Vec<f64> = vec![0.0; n];
        let mut gbest_x: Genome = Genome::new(n + 1, d);
        let mut gbesthistory: Vec<f64> = vec![0.0; max_iter];

        let mut best_x: Genome = Genome::new(n + 2, d);
        let mut worst_x: Genome = Genome::new(n + 3, d);
        let mut better_x: Genome = Genome::new(n + 4, d);
        let mut r: Genome = Genome::new(n + 5, d);

        let mut gap: Vec<Vec<f64>> = vec![vec![0.0; d]; 4];
        let mut distance: [f64; 4] = [0.0; 4];
        let mut lf: [f64; 4] = [0.0; 4];
        let mut ka: Vec<Vec<f64>> = vec![vec![0.0; d]; 4];

        let mut rng = rand::thread_rng();
        let intervall01 = Uniform::from(0.0..=1.0);
        let intervall0_p1 = Uniform::from(0..P1);

        let mut new_x: Vec<Genome> = self.get_empty_solutions(n);

        //Initialization
        let mut x = self.initialize(self.params, InitializationMode::RealUniform);

        //Evaluation of search agents
        // Sequential mode
        #[cfg(not(feature = "parallel"))]
        for i in 0..n {
            fitness[i] = self.problem.objectivefunction(&x[i].genes);
            fes += 1.0;
        }

        //___________Parallel mode________________
        #[cfg(feature = "parallel")]
        {
            x.par_iter_mut()
                .for_each(|g| g.fitness = Some(self.problem.objectivefunction(&g.genes)));
            for i in 0..n {
                match x[i].fitness {
                    None => fitness[i] = f64::MAX,
                    Some(fit) => fitness[i] = fit,
                };
            }
        }
        //________________________________________

        // Save the best solution
        for i in 0..n {
            if gbestfitness > fitness[i] {
                gbestfitness = fitness[i];
                copy_vector(&x[i].genes, &mut gbest_x.genes);
            }
        }

        gbesthistory[0] = gbestfitness;

        //println!("Best_fitness : {}", gbestfitness);

        while iter < max_iter {
            //Sorte and sorting indexes:
            let mut ind: Vec<usize> = (0..n).collect();
            ind.sort_by(|&a, &b| fitness[a].total_cmp(&fitness[b]));
            //----------------------------------------------------------------------------------------

            //println!("ind : {:?}", ind);

            // Save best solution
            copy_vector(&x[ind[0]].genes, &mut best_x.genes);

            // Learning phase
            let interval_worst = Uniform::from(n - P1..n);
            let interval_better = Uniform::from(1..P1);

            for i in 0..n {
                //  Worst_X = x(ind(randi([popsize-P1+1,popsize],1)),:);
                let wors_index: usize = interval_worst.sample(&mut rng);
                copy_vector(&x[ind[wors_index]].genes, &mut worst_x.genes);

                // Better_X=x(ind(randi([2,P1],1)),:);
                let better_index: usize = interval_better.sample(&mut rng);
                copy_vector(&x[ind[better_index]].genes, &mut better_x.genes);

                //random=selectID(popsize,i,2); L1=random(1); L2=random(2);
                let (l1, l2) = self.select_id(i, n, &mut rng);

                // Gap1=(Best_X-Better_X);
                for j in 0..d {
                    gap[0][j] = best_x.genes[j] - better_x.genes[j];
                }

                //Gap2=(Best_X-Worst_X);
                for j in 0..d {
                    gap[1][j] = best_x.genes[j] - worst_x.genes[j];
                }

                //Gap3=(Better_X-Worst_X);
                for j in 0..d {
                    gap[2][j] = better_x.genes[j] - worst_x.genes[j];
                }

                //Gap4=(x(L1,:)-x(L2,:));
                for j in 0..d {
                    gap[3][j] = x[l1].genes[j] - x[l2].genes[j];
                }

                //Distance1=norm(Gap1); Distance2=norm(Gap2); Distance3=norm(Gap3); Distance4=norm(Gap4);
                for k in 0..4 {
                    distance[k] = self.norm(&gap[k]);
                }

                // SumDistance=Distance1+Distance2+Distance3+Distance4;
                let mut sum_distance = distance.iter().fold(0.0, |sum, a| sum + a);
                if sum_distance == 0.0 {
                    sum_distance = 1.0;
                }

                // LF1=Distance1/SumDistance;  LF2=Distance2/SumDistance;  LF3=Distance3/SumDistance; LF4=Distance4/SumDistance;

                for k in 0..4 {
                    lf[k] = distance[k] / sum_distance;
                }

                let max_fitness = match fitness.iter().max_by(|a, b| a.total_cmp(b)) {
                    Some(value) => *value,
                    None => fitness[0],
                };

                //  SF=(fitness(i)/max(fitness));
                let sf = fitness[i] / max_fitness;

                //KA1=LF1*SF*Gap1; KA2=LF2*SF*Gap2; KA3=LF3*SF*Gap3; KA4=LF4*SF*Gap4;
                for t in 0..4 {
                    for j in 0..d {
                        ka[t][j] = sf * lf[t] * gap[t][j];
                    }
                }

                // newx(i,:)=x(i,:)+KA1+KA2+KA3+KA4;
                //let sum_ka = ka.iter().fold(0.0, |sum, a|  sum + a);
                for j in 0..d {
                    new_x[i].genes[j] = x[i].genes[j] + ka[0][j] + ka[1][j] + ka[2][j] + ka[3][j];
                }

                // Space bound
                for j in 0..d {
                    new_x[i].genes[j] = f64::min(new_x[i].genes[j], ub[j]);
                    new_x[i].genes[j] = f64::max(new_x[i].genes[j], lb[j]);
                }

                //newfitness= ObjectiveFunction(newx(i,:));

                let new_fitness = self.problem.objectivefunction(&new_x[i].genes);
                fes += 1.0;

                //Update solutions
                if new_fitness < fitness[i] {
                    fitness[i] = new_fitness;
                    copy_solution(&new_x[i], &mut x[i], d);
                } else {
                    let rand_value = intervall01.sample(&mut rng);
                    if rand_value < P2 && ind[i] != ind[0] {
                        fitness[i] = new_fitness;
                        copy_solution(&new_x[i], &mut x[i], d);
                    }
                }

                // Save the best solution
                if gbestfitness > fitness[i] {
                    gbestfitness = fitness[i];
                    copy_solution(&x[i], &mut gbest_x, d);
                }
            }

            // Reflection phase
            for i in 0..n {
                //newx(i,:)=x(i,:);
                copy_solution(&x[i], &mut new_x[i], d);

                for j in 0..d {
                    let rand_value = intervall01.sample(&mut rng);
                    if rand_value < P3 {
                        // R=x(ind(randi(P1)),:);
                        let rand_index = intervall0_p1.sample(&mut rng);
                        copy_solution(&x[ind[rand_index]], &mut r, d);

                        // newx(i,j) = x(i,j)+(R(:,j)-x(i,j))*unifrnd(0,1);
                        new_x[i].genes[j] = x[i].genes[j] + r.genes[j]
                            - x[i].genes[j] * intervall01.sample(&mut rng);

                        // AF=(0.01+(0.1-0.01)*(1-FEs/MaxFEs));
                        let af = 0.01 + (0.1 - 0.01) * (1.0 - fes / max_fes);

                        if intervall01.sample(&mut rng) < af {
                            //newx(i,j)=xmin+(xmax-xmin)*unifrnd(0,1);
                            new_x[i].genes[j] =
                                lb[j] + (ub[j] - lb[j]) * intervall01.sample(&mut rng);
                        }
                    }
                }

                // Space bound
                for j in 0..d {
                    new_x[i].genes[j] = f64::min(new_x[i].genes[j], ub[j]);
                    new_x[i].genes[j] = f64::max(new_x[i].genes[j], lb[j]);
                }
                //___________________________________________________________________________
                //  newfitness= ObjectiveFunction(newx(i,:));
                let new_fitness = self.problem.objectivefunction(&new_x[i].genes);
                //FEs=FEs+1;
                fes += 1.0;

                //Update solutions
                if new_fitness < fitness[i] {
                    fitness[i] = new_fitness;
                    copy_solution(&new_x[i], &mut x[i], d);
                } else {
                    let rand_value = intervall01.sample(&mut rng);
                    if rand_value < P2 && ind[i] != ind[0] {
                        fitness[i] = new_fitness;
                        copy_solution(&new_x[i], &mut x[i], d);
                    }
                }
                // Save the best solution
                if gbestfitness > fitness[i] {
                    gbestfitness = fitness[i];
                    copy_solution(&x[i], &mut gbest_x, d);
                }
            }

            gbesthistory[iter] = gbestfitness;
            #[cfg(feature = "report")]
            println!("Iter : {}, best-fitness : {}", iter, gbestfitness);

            iter += 1;
        }

        let duration = chronos.elapsed();
        //println!("Iter : {}, FES: {}, Max_FES {}", iter, fes, max_fes);

        // Compute the best fitness for the best solution
        best_x.fitness = Some(self.problem.objectivefunction(&best_x.genes));

        let result: OptimizationResult = OptimizationResult {
            best_genome: Some(best_x),
            best_fitness: Some(gbestfitness),
            convergence_trend: Some(gbesthistory),
            computation_time: Some(duration),
            err_report: None,
        };

        result
    }
}

///
/// Define parameters for the Growth (GO) algorithm.
///
#[derive(Debug, Clone)]
pub struct GOparams<'a> {
    /// The number of search agents.
    pub population_size: usize,

    /// The dimension of the optimization problem (i.e., the length of the solution).
    pub problem_dimension: usize,

    /// The maximum number of iterations serves as the stopping criterion for the optimization process.
    pub max_iterations: usize,

    /// The lower bounds of the search space.
    pub lower_bounds: &'a [f64],

    /// The upper bounds of the search space.
    pub upper_bounds: &'a [f64],
}

impl<'a> GOparams<'a> {
    ///
    /// Create a new instance of GO parameters:
    /// pop_size : The number of search agents.
    /// dim : The dimension of the optimization problem.
    /// max_iter : The maximum number of iterations serves as the stopping criterion.
    /// lb : The lower bounds of the search space.
    /// ub : The upper bounds of the search space.
    ///
    #[allow(dead_code)]
    pub fn new(pop_size: usize, dim: usize, max_iter: usize, lb: &'a [f64], ub: &'a [f64]) -> Self {
        GOparams {
            population_size: pop_size,
            problem_dimension: dim,
            max_iterations: max_iter,
            lower_bounds: lb,
            upper_bounds: ub,
        }
    }
}

impl<'a> Parameters for GOparams<'a> {
    fn get_problem_dimension(&self) -> usize {
        self.problem_dimension
    }

    fn get_max_iterations(&self) -> usize {
        usize::max(self.max_iterations, 1)
    }

    fn get_population_size(&self) -> usize {
        usize::max(self.population_size, 6)
    }

    fn get_lower_bounds(&self) -> Vec<f64> {
        self.lower_bounds.to_vec()
    }

    fn get_upper_bounds(&self) -> Vec<f64> {
        self.upper_bounds.to_vec()
    }
}

impl<'a> Default for GOparams<'a> {
    ///
    /// Return the default values of parameters, as follows:
    ///
    /// ~~~
    ///
    ///  use sefar::algos::go::*;
    ///
    ///  GOparams {
    ///     population_size : 20,
    ///     dimensions : 3,
    ///     max_iterations : 100,
    ///     lower_bounds : &[100.0f64, 100.0, 100.0],
    ///     upper_bounds : &[-100.0f64, -100.0, -100.0],
    /// };
    /// ~~~
    ///
    fn default() -> Self {
        GOparams {
            population_size: 20,
            problem_dimension: 3,
            max_iterations: 100,
            lower_bounds: &[-100.0f64, -100.0, -100.0],
            upper_bounds: &[100.0f64, 100.0, 100.0],
        }
    }
}

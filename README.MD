# Sefar 

[Sefar](https://github.com/SaadDAHMANI/sefar) is a simple and comprehensive [Rust](https://github.com/rust-lang/rust) library for evolutionary optimization algorithms, exclusively written using safe code. It supports **continuous** and **binary** optimization in both **sequential** and **parallel** modes through its features. In the current version, the *_parallel mode executes objective function_* evaluations in parallel (multi-threading) using the [rayon](https://github.com/rayon-rs/rayon) crate.

## Current state (Under development)

 In this version, [Sefar](https://github.com/SaadDAHMANI/sefar) supports: 

- [X] Particle Swarm Optimization ([PSO](https://doi.org/10.1109/ICNN.1995.488968));   
- [X] Equilibrium optimizer ([EO](https://doi.org/10.1016/j.knosys.2019.105190));
- [X] Growth Optimizer ([GO](https://doi.org/10.1016/j.knosys.2022.110206)).

## Binary optimization
In the current version, binarization is performed using the S-Shape function provided below:

$S(x) = 1/(1 + e^{(-x)})$

In this context, *x* represents a "gene" and signifies each element in the candidate solution *X* ("genome") within a search space of length *d*, where $X= \{x_1, x_2, ..., x_d\}$. 

The Binary optimization can be executed using the **binary** feature.

### Example

```Rust
use sefar::core::eoa::EOA;
use sefar::core::optimization_result::OptimizationResult;
use sefar::algos::go::{GOparams, GO};

fn main() {
    println!("Hello, sefar !");
   
    #[cfg(feature ="binary")] go_f1_binary_test();
}

///
/// run the binary version of Growth Optimizer (Binary-GO).
/// 
#[cfg(feature = "binary")]
fn go_f1_binary_test(){

    // Define the parameters of GO:
    let search_agents : usize = 20;
    let dim : usize = 10;
    let max_iterations : usize = 50;
    let lb = vec![0.0; dim];
    let ub = vec![1.0; dim];
    
    // Build the parameter struct:
    let mut settings : GOparams = GOparams::new(search_agents, dim, max_iterations, &lb, &ub);
    
    // Define the problem to optimize:
    let mut fo = Sphere{};
  
    // Build the optimizer:
    let mut algo : GO<Sphere> = GO::new(&settings, &mut fo);
    
    // Run the GO algorithm: 
    let result : OptimizationResult = algo.run();

    // Print the results:
    println!("The optimization results of Binary-GO : {}", result.to_string());
} 

///
/// Sphere benchmark function (F1). 
/// Fi(X) = Sum(|X|)
/// where X = {x1, x2, ..... xd}, and 'd' is the problem dimension.
/// 
#[derive(Debug,Clone)]
pub struct Sphere{}

impl Problem for Sphere{
    #[cfg(not(feature="parallel"))]
    fn objectivefunction(&mut self, genome : &[f64]) ->f64 {
        let fitness = genome.iter().fold(0.0f64, |sum, g| sum + g.powi(2));
        fitness        
    }

    #[cfg(feature="parallel")]
    fn objectivefunction(&self, genome : &[f64]) ->f64 {
        let fitness = genome.iter().fold(0.0f64, |sum, g| sum + g.powi(2));
        fitness        
    }   
}

```


## Supported features

|Feature        | Designation                                         |
| ------------- | --------------------------------------------------- |
| *_binary_*    | run binary optimization using **S-Shape** function  |
| *_parallel_*  | run optimization in parallel mode using Rayon crate | 


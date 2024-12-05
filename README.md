# Sefar

[Sefar](https://github.com/SaadDAHMANI/sefar) is a simple and comprehensive [Rust](https://github.com/rust-lang/rust) library for evolutionary optimization algorithms, exclusively written using Rust safe code. It supports **continuous** and **binary** optimization in both **sequential** and **parallel** modes. In the current version, the *_parallel mode executes objective function_* evaluations in parallel (multi-threading) using [rayon](https://github.com/rayon-rs/rayon) crate.

## Current state (Under development)

1. [Sefar](https://github.com/SaadDAHMANI/sefar) perfoms **minimization** by default. In the case of **maximization**, the objective function $f(X)$ can be expressed as $-f(X)$.

2. In this version, [Sefar](https://github.com/SaadDAHMANI/sefar) supports:

- [X] Particle Swarm Optimization ([PSO](https://doi.org/10.1109/ICNN.1995.488968));
- [X] Equilibrium optimizer ([EO](https://doi.org/10.1016/j.knosys.2019.105190));
- [X] Binary Equilibrium Optimizer ([BiEO](https://doi.org/10.1016/j.enbuild.2022.112503));
- [-] Modified Equilibrium optimizer ([MEO](https://doi.org/10.1016/j.asoc.2020.106542));
- [X] Growth Optimizer ([GO](https://doi.org/10.1016/j.knosys.2022.110206));
- [X] Gaining-Sharing Knowledge ([GSK](https://doi.org/10.1007/s13042-019-01053-x));
- [-] Gaining-Sharing Knowledge Based Algorithm With Adaptive Parameters ([APGSK](https://doi:10.1109/ACCESS.2021.3076091));
- [-] LSHADE_SPACMA([LSHADE_SPACMA](https://ieeexplore.ieee.org/document/7969307))

## Binary optimization
The binary optimizatin in the older versions of [Sefar](https://github.com/SaadDAHMANI/sefar) will be replaced by more efficent binary optimization algorithms.*
This approach aims to simplify the implementation of algorithms on one hand and to offer parallel mode for binary algorithms as well in simple way.

## Parallel optimization
In the current version of [Sefar](https://github.com/SaadDAHMANI/sefar), only the objective function evaluation is run in parallel mode using [rayon](https://github.com/rayon-rs/rayon) crate.

### Example
1. Import [Sefar](https://github.com/SaadDAHMANI/sefar) in the *Cargo.Toml* file of your project.

```toml

[dependencies]
sefar = "0.1.7"
```

2. In the *main.rs* file :

```rust
extern crate sefar;
use sefar::core::eoa::EOA;
use sefar::core::optimization_result::OptimizationResult;
use sefar::algos::go::{GOparams, GO};
use sefar::core::problem::Problem;

fn main() {

    println!("Optimization using Growth optimizer in Sefar crate:");

   // Define the parameters of GO:
    let search_agents : usize = 20; // number of search agents.
    let dim : usize = 5; // problem dimension.
    let max_iterations : usize = 200; // maximum number of iterations.
    let lb = vec![-100.0; dim]; // lower bound of search space.
    let ub = vec![100.0; dim]; // upper bound of the search space.

    // Build the parameter struct:
    let settings : GOparams = GOparams::new(search_agents, dim, max_iterations, &lb, &ub);

    // Define the problem to optimize:
    let mut fo = F1{};

    // Build the optimizer:
    let mut algo : GO<F1> = GO::new(&settings, &mut fo);

    // Run the GO algorithm:
    let result : OptimizationResult = algo.run();

    // Print the results:
    println!("The optimization results of GO : {}", result.to_string());

    // The result will be something like :
    // The optimization results of GO : Best-fitness : Some(1.2106003206412792e-54);
    // Best-solution : Some(Genome { id: 22, genes: [-7.586125521377413e-28, -7.519595439155215e-28, -2.2218253597758204e-29, -6.135485510888784e-29, -3.7827445210037567e-28], fitness: Some(1.2882857900967827e-54) });
    // Time : Some(11.606888ms);
    // Err-report: None
}

// Define the objective function to minimize. Here, the Sphere function is implemented.

///
/// F1 : Sphere benchmark function.
/// Fi(X) = Sum(|X^2|)
/// where X = {x1, x2, ..... xd}, and 'd' is the problem dimension.
///
#[derive(Debug,Clone)]
pub struct F1{}

impl Problem for F1 {
    fn objectivefunction(&mut self, genome : &[f64])->f64 {
       genome.iter().fold(0.0f64, |sum, x| sum + x.powi(2))
    }
}
```

## Supported features
[Sefar](https://github.com/SaadDAHMANI/sefar) supports *_report_* and *_parallel_* features.

1. Run *_report_* feature:

```bash
# run report feature:
cargo run --features report;

# run parallel feature:
cargo run --features parallel;
```

|Algorithm       | *_report_* feature | *_parallel_* feature |
|----------------|--------------------| ---------------------|
|*PSO*           | :heavy_check_mark: |                      |
|*EO*            | :heavy_check_mark: | :heavy_check_mark:   |
|*GO*            | :heavy_check_mark: | :heavy_check_mark:   |
|*GSK*           | :heavy_check_mark: | :heavy_check_mark:   |
|*LSHADE_SPACMA* |                    |                      |
|*BiEO*          | :heavy_check_mark: |                      |

# Sefar

[Sefar](https://github.com/SaadDAHMANI/sefar) is a simple and comprehensive [Rust](https://github.com/rust-lang/rust) library for evolutionary optimization algorithms, exclusively written using Rust safe code. It supports **continuous** and **binary** optimization in both **sequential** and **parallel** modes through its features. In the current version, the *_parallel mode executes objective function_* evaluations in parallel (multi-threading) using the [rayon](https://github.com/rayon-rs/rayon) crate.

## Current state (Under development)

1. [Sefar](https://github.com/SaadDAHMANI/sefar) perfoms **minimization** by default. In the case of **maximization**, the objective function $f(X)$ can be expressed as $-f(X)$.

2. In this version, [Sefar](https://github.com/SaadDAHMANI/sefar) supports:

- [X] Particle Swarm Optimization ([PSO](https://doi.org/10.1109/ICNN.1995.488968));
- [X] Equilibrium optimizer ([EO](https://doi.org/10.1016/j.knosys.2019.105190));
- [-] Modified Equilibrium optimizer ([MEO](https://doi.org/10.1016/j.asoc.2020.106542));
- [X] Growth Optimizer ([GO](https://doi.org/10.1016/j.knosys.2022.110206));
- [X] Gaining-Sharing Knowledge ([GSK](https://doi.org/10.1007/s13042-019-01053-x));
- [-] Gaining-Sharing Knowledge Based Algorithm With Adaptive Parameters ([APGSK](https://doi:10.1109/ACCESS.2021.3076091));
- [-] LSHADE_SPACMA([LSHADE_SPACMA](https://ieeexplore.ieee.org/document/7969307))


# Important

**In the current version, binary and parallel optimization are implemented exclusively for the Equilibrium Optimizer (EO) and the Growth Optimizer (GO). Soon, these features will be available for the other algorithms as well.**

## Binary optimization

In the current version, binarization is performed using the S-Shape function provided below:

$S(x) = 1/(1 + e^{(-x)})$

In this context, *x* represents a "gene" and signifies each element in the candidate solution *X* ("genome") within a search space of length *d*, where $X= \{x_1, x_2, ..., x_d\}$.

The Binary optimization can be executed using the **binary** feature.

### Example
1. Import [Sefar](https://github.com/SaadDAHMANI/sefar) with **binary** feature in the *Cargo.Toml* file of your project.

```toml

[dependencies]
sefar = {version = "0.1.3", features = ["binary"]}
```

2. In the *main.rs* file :

```rust
extern crate sefar;
use sefar::core::eoa::EOA;
use sefar::core::optimization_result::OptimizationResult;
use sefar::algos::go::{GOparams, GO};
use sefar::core::problem::Problem;

fn main() {

    println!("Binary optimization using Growth optimizer in Sefar crate:");

    go_f1_binary_test();
}

///
/// run the binary version of Growth Optimizer (Binary-GO).
///
fn go_f1_binary_test(){

    // Define the parameters of GO:
    let search_agents : usize = 20;
    let dim : usize = 10;
    let max_iterations : usize = 100;
    let lb = vec![0.0; dim];
    let ub = vec![1.0; dim];

    // Build the parameter struct:
    let settings : GOparams = GOparams::new(search_agents, dim, max_iterations, &lb, &ub);

    // Define the problem to optimize:
    let mut fo = F1{};

    // Build the optimizer:
    let mut algo : GO<F1> = GO::new(&settings, &mut fo);

    // Run the GO algorithm:
    let result : OptimizationResult = algo.run();

    // Print the results:
    println!("The optimization results of Binary-GO : {}", result.to_string());

    // The results show something like :
    // Binary optimization using Growth optimizer in Sefar crate:
    // The optimization results of Binary-GO : Best-fitness : Some(0.0)
    // ; Best-solution : Some(Genome { id: 22, genes: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], fitness: Some(0.0) })
    // ; Time : Some(3.326498ms)
    // ; Err-report: None
}

// Define the objective function to minimize. Here, the Sphere function is implemented.

///
/// F1 : Sphere benchmark function.
/// Fi(X) = Sum(|X|)
/// where X = {x1, x2, ..... xd}, and 'd' is the problem dimension.
///
#[derive(Debug,Clone)]
pub struct F1{}

impl Problem for F1 {
    fn objectivefunction(&mut self, genome : &[f64])->f64 {
       genome.iter().fold(0.0f64, |sum, x| sum +x)
    }
}
```

## Supported features

To run *report* feature:

```bash
cargo run --features report;
```

|Algorithm       | *_report_* feature | *_binary_* feature   |  *_parallel_* feature |
|*PSO*           | [x]                | [ ]                  | [ ]                   |
|*EO*            | [x]                | [x] S-Shape function | [x]                   |
|*GO*            | [x]                | [x] S-Shape functio  | [x]                   |
|*GSK*           | [x]                | [ ]                  | [ ]                   |
|*LSHADE_SPACMA* | [x]                | [ ]                  | [ ]                   |

# Sefar 

[Sefar](https://github.com/SaadDAHMANI/sefar) is a simple and comprehensive [Rust](https://github.com/rust-lang/rust) library for evolutionary optimization algorithms, exclusively written using safe code. It supports **continuous** and **binary** optimization in both **sequential** and **parallel** modes through its features. In the current version, the *_parallel mode executes objective function_* evaluations in parallel (multi-threading) using the [rayon](https://github.com/rayon-rs/rayon) crate.

## Current state (Under development)

1. [Sefar](https://github.com/SaadDAHMANI/sefar) perfoms **minimization** by default. In the case of **maximization**, the objective function $f(X)$ can be expressed as $-f(X)$.

2. In this version, [Sefar](https://github.com/SaadDAHMANI/sefar) supports: 

- [X] Particle Swarm Optimization ([PSO](https://doi.org/10.1109/ICNN.1995.488968));   
- [X] Equilibrium optimizer ([EO](https://doi.org/10.1016/j.knosys.2019.105190));
- [X] Growth Optimizer ([GO](https://doi.org/10.1016/j.knosys.2022.110206)).

## Binary optimization
In the current version, binarization is performed using the S-Shape function provided below:

$S(x) = 1/(1 + e^{(-x)})$

In this context, *x* represents a "gene" and signifies each element in the candidate solution *X* ("genome") within a search space of length *d*, where $X= \{x_1, x_2, ..., x_d\}$. 

The Binary optimization can be executed using the **binary** feature.

### Example
1. Import [Sefar](https://github.com/SaadDAHMANI/sefar) with **binary** feature in the *Cargo.Toml* file of your project.

2. In the *main.rs* file :

```Rust


```


## Supported features

|Feature        | Designation                                         |
| ------------- | --------------------------------------------------- |
| *_binary_*    | run binary optimization using **S-Shape** function  |
| *_parallel_*  | run optimization in parallel mode using Rayon crate | 


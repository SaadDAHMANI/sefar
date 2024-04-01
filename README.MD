# Sefar 

[Sefar](https://github.com/SaadDAHMANI/sefar) is a simple and comprehensive [Rust](https://github.com/rust-lang/rust) library for evolutionary optimization algorithms, exclusively written using safe code. It supports **continuous** and **binary** optimization in both **sequential** and **parallel** modes through its features. In the current version, the *_parallel mode executes objective function_* evaluations in parallel (multi-threading) using the [rayon](https://github.com/rayon-rs/rayon) crate.

## Current state (Under development)

 Sefar supports currently : 

- [X] Particle Swarm Optimization ([PSO](https://doi.org/10.1109/ICNN.1995.488968));   
- [X] Equilibrium optimizer ([EO](https://doi.org/10.1016/j.knosys.2019.105190));
- [X] Growth Optimizer ([GO](https://doi.org/10.1016/j.knosys.2022.110206)).

## Binary optimization
In the current version 


## Supported features

|Feature        | Designation  |
| ------------- | --------------------- |
| *_binary_*    | run binary optimization using **S-Shape** function. Use : *_cargo run --features "binary"_* | 
| *_parallel_*  | run optimization in parallel mode | 


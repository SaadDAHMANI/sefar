//
// Implementation of Equilibrium Optimizer
// 
 
extern crate rand;
use rand::distributions::{Distribution, Uniform};
use std::time::Instant;

use crate::core::genome::Genome;

use crate::core::problem::Problem;
use crate::core::optimization_result::OptimizationResult;
//use crate::sequential_algos::common::*;
use crate::common::*;


pub struct EO<'a, T : Problem> {
    problem : &'a mut T,
}
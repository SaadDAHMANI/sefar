//
// Implementation of Particle Swarm Optimization algorithm (PSO)
// 
 
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

pub struct PSO<'a, T: Problem> {
    problem : &'a mut T,
     
}
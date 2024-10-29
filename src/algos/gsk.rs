use crate::common::*;
use crate::core::eoa::{InitializationMode, EOA};
use crate::core::genome::Genome;
use crate::core::optimization_result::OptimizationResult;
use crate::core::parameters::Parameters;
use crate::core::problem::Problem;
use rand::rngs::ThreadRng;
use rand_distr::{Distribution, Uniform};
//#[cfg(feature = "parallel")]
//use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use std::time::Instant;

/// GSK : Gaining-Sharing Knowedge algorithm.
/// Reference:
/// Mohamed, A. W., Hadi, A. A., & Mohamed, A. K. (2020).
/// Gaining-sharing knowledge based algorithm for solving optimization problems: a novel nature-inspired algorithm.
/// International Journal of Machine Learning and Cybernetics, 11(7), 1501-1529.
/// (https​://doi.org/10.1007/s1304​2-019-01053​-x)
///

pub mod eoa;
pub mod genome;
pub mod optimization_result;
pub mod parameters;
pub mod problem;

use optimization_result::OptimizationResult;
use thiserror::*;

pub trait EvolutionaryAlgo {
    fn run_epoch(&mut self) {}
    fn run(&mut self) -> Option<OptimizationResult> {
        None
    }
}

///
/// A custom error to represent EOA errors.
///
#[derive(Debug, Clone, Error, PartialEq)]
pub enum OptError {
    #[error("Search population size is eqal to 0!")]
    PopulationSizeIsNull,

    #[error("Search population size is less than the minimum required!")]
    PopulationSizeLessMin { min: usize, actual: usize },

    #[error("Problem dimension is equal to 0!")]
    ProblemDimensionIsNull,

    #[error("Maximum number of iterations is equal to 0!")]
    MaxIterationsIsNull,

    #[error("Length of Lower-bound vector is not equal to problem dimension!")]
    LBLengthNoEqualsProblemDim { lb_len: usize, problem_dim: usize },

    #[error("Length of upper-bound vector is not equal to problem dimension!")]
    UBLengthNoEqualsProblemDim { ub_len: usize, problem_dim: usize },

    #[error("Lengths of Lower-bound and upper-bound vectors are not equal!")]
    LBLengthNotEqualsUBlength { lb_len: usize, ub_len: usize },

    #[error("Lower-bound vector is empty!")]
    EmptyLB,

    #[error("Upper-bound vector is empty!")]
    EmptyUB,

    #[error("Bad parameter value.")]
    BadParameterValue { parameter: String, actual: f64 },

    #[error(" ")]
    ThreadPoolBuildErr,

    #[error("Undefined error")]
    Other,
}

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
#[derive(Debug, Clone, Copy, Error, PartialEq, Eq)]
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
    LBLengthNoEqualsProblemDim,

    #[error("Length of upper-bound vector is not equal to problem dimension!")]
    UBLengthNoEqualsProblemDim,

    #[error("Lengths of Lower-bound and upper-bound vectors are not equal!")]
    LBLengthNotEqualsUBlength,

    #[error("Lower-bound vector is empty!")]
    EmptyLB,

    #[error("Upper-bound vector is empty!")]
    EmptyUB,
}

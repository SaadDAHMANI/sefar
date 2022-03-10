pub mod benchmarks;
pub mod core;
pub mod sequential_algos;
//mod parallel_algos;
mod common;
//mod paracommon;

use crate::benchmarks::functions::*;
use crate::sequential_algos::eo::*;
//use crate::sequential_algos::pso::*;
//use crate::sequential_algos::ga::*;



#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}

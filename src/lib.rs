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
    use crate::sequential_algos::eo::EOparams;

    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }


    #[test]
    fn eo_f1_test1(){
        let settings : EOparams = EOparams::default();
        
    }
}

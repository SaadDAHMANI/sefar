pub mod benchmarks;
pub mod core;
pub mod sequential_algos;
//mod parallel_algos;
mod common;
//mod paracommon;


#[cfg(test)]
mod tests {
    use crate::benchmarks::functions::F1;
    use crate::sequential_algos::eo::EO;
    use crate::sequential_algos::eo::EOparams;
    use crate::core::eoa::EOA;
    //use super::*;

    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }


    #[test]
    fn eo_f1_test1(){
        let settings : EOparams = EOparams::default();
        //settings.population_size = 20;
        //settings.max_iterations = 500; 

        let mut fo = F1 {};

        let mut eo : EO<F1> = EO::new(&settings, &mut fo);
        
        let result = eo.run();

        let bestfit = result.best_fitness.unwrap();


        assert_eq!(result.best_fitness, Some(0.0001));
        
     }
}

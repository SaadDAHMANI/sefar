pub mod algos;
pub mod benchmarks;
pub mod core;
//pub mod parallel;
//mod parallel_algos;
mod common;

//mod paracommon;

#[cfg(test)]
mod tests {
    use crate::algos::eo::eo::EOparams;
    use crate::algos::eo::eo::EO;
    use crate::algos::pso::*;
    use crate::benchmarks::functions::Sphere;
    use crate::core::eoa::EOA;
    //use super::*;

    #[test]
    fn eo_f1_test1() {
        let settings: EOparams = EOparams::default();
        //settings.population_size = 20;
        //settings.max_iterations = 500;

        let mut fo = Sphere {};

        let mut eo: EO<Sphere> = EO::new(&settings, &mut fo);

        let result = eo.run();

        //let bestfit = result.best_fitness.unwrap();

        assert_ne!(result.best_fitness, None);
    }

    #[test]
    fn pso_f1_test1() {
        let settings: PSOparams = PSOparams::default();
        //settings.population_size = 20;
        //settings.max_iterations = 500;

        let mut fo = Sphere {};

        let mut pso: PSO<Sphere> = PSO::new(&settings, &mut fo);

        let result = pso.run();

        //let bestfit = result.best_fitness.unwrap();

        assert_ne!(result.best_fitness, None);
    }
}

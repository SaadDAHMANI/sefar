//use crate::core::problem::Objectivefunction;
use std::cmp::Ordering;

#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub struct Genome {
    pub id: usize,
    pub genes: Vec<f64>,
    pub fitness: Option<f64>,
}
#[allow(dead_code)]
impl Genome {
    //
    // Retun new Genome wit fitness = f64::MAX
    //
    //
    pub fn new(genome_id: usize, dimension: usize) -> Genome {
        Genome {
            id: genome_id,
            genes: vec![0.0f64; dimension],
            fitness: Some(f64::MAX),
        }
    }

    pub fn get_dimensions(&self) -> usize {
        self.genes.len()
    }

    pub fn from(the_id: usize, solution: &[f64], fitnes: f64) -> Genome {
        Genome {
            id: the_id,
            genes: solution.to_owned(),
            fitness: Some(fitnes),
        }
    }

    pub fn cmp_genome(a: &Genome, b: &Genome) -> Ordering {
        if a.fitness == None && b.fitness == None {
            return Ordering::Equal;
        }

        if a.fitness == None && b.fitness != None {
            return Ordering::Greater;
        }

        if a.fitness != None && b.fitness == None {
            return Ordering::Less;
        }

        if a.fitness.unwrap().is_nan() {
            return Ordering::Greater;
        }

        if b.fitness.unwrap().is_nan() {
            return Ordering::Less;
        }

        if a.fitness.unwrap() < b.fitness.unwrap() {
            return Ordering::Less;
        } else if a.fitness.unwrap() > b.fitness.unwrap() {
            return Ordering::Greater;
        }

        Ordering::Equal
    }
}

#[cfg(test)]
mod genome_tests {
    use super::*;

    #[test]
    fn cmp_genome_test1() {
        let mut gen1 = Genome::new(0, 3);
        let mut gen2 = Genome::new(1, 3);
        let mut gen3 = Genome::new(2, 3);
        let mut gen4 = Genome::new(3, 3);

        gen1.fitness = Some(0.5f64);
        gen2.fitness = Some(0.0001f64);
        gen3.fitness = Some(50.5f64);
        gen4.fitness = None;

        let mut pop = Vec::new();
        pop.push(gen1);
        pop.push(gen2);
        pop.push(gen3);
        pop.push(gen4);

        pop.sort_by(Genome::cmp_genome);

        let mut correct = Vec::new();
        correct.push(pop[1].clone());
        correct.push(pop[0].clone());
        correct.push(pop[2].clone());
        correct.push(pop[3].clone());

        assert_ne!(correct, pop);
    }
}

extern crate rand;
extern crate rand_distr;
//#[cfg(feature = "binary")]
//use rand::rngs::ThreadRng;

use rand::distributions::{Distribution, Uniform};
use rand_distr::Normal;

use crate::core::genome::Genome;

#[allow(dead_code)]
pub fn copy_matrix(source: &[Genome], destination: &mut Vec<Vec<f64>>) {
    let ni = source.len();
    let nj = source[0].get_dimensions();

    for i in 0..ni {
        for j in 0..nj {
            destination[i][j] = source[i].genes[j];
        }
    }
}

pub fn randomize(randvect: &mut Vec<f64>) {
    let between = Uniform::from(0.0..=1.0);
    let mut rng = rand::thread_rng();

    for i in 0..randvect.len() {
        randvect[i] = between.sample(&mut rng);
    }
}

#[allow(dead_code)]
pub fn uniform_rand(min_value: f64, max_value: f64) -> Result<f64, String> {
    if min_value > max_value {
        Err(String::from("uniform_rand : min_value > max_value !!"))
    } else {
        let between = Uniform::from(min_value..=max_value);
        let mut rng = rand::thread_rng();
        Ok(between.sample(&mut rng))
    }
}

#[allow(dead_code)]
pub fn uniform_rand1(min_value: f64, max_value: f64) -> f64 {
    let between = Uniform::from(min_value..=max_value);
    let mut rng = rand::thread_rng();
    between.sample(&mut rng)
}

#[allow(dead_code)]
pub fn normal_rand1(mean: f64, std_dev: f64) -> f64 {
    let normal = Normal::new(mean, std_dev).unwrap();
    let v = normal.sample(&mut rand::thread_rng());
    v
}

/// Return a vector with Uniform random distribution
#[allow(dead_code)]
pub fn rand_vec(imin: usize, imax: usize, cols: usize) -> Vec<usize> {
    let mut vector: Vec<usize> = vec![0; cols];

    //println!("i_min {}, imax : {}", imin, imax);
    if imin == imax {
        for j in 0..cols {
            vector[j] = imin;
        }
    } else {
        let min_u = usize::min(imin, imax);
        let max_u = usize::max(imin, imax);

        let between = Uniform::from(min_u..=max_u);
        let mut rng = rand::thread_rng();

        for j in 0..cols {
            vector[j] = between.sample(&mut rng);
        }
    }
    vector
}

///
/// Return a vector with Uniform random distribution
///
#[allow(dead_code)]
pub fn rand_matrix(min_value: f64, max_value: f64, rows: usize, cols: usize) -> Vec<Vec<f64>> {
    let mut matrix: Vec<Vec<f64>> = vec![vec![0.0; cols]; rows];

    if min_value >= max_value {
        for i in 0..rows {
            for j in 0..cols {
                matrix[i][j] = max_value;
            }
        }
    } else {
        let between = Uniform::from(min_value..max_value);
        let mut rng = rand::thread_rng();
        for i in 0..rows {
            for j in 0..cols {
                matrix[i][j] = between.sample(&mut rng);
            }
        }
    }
    matrix
}

pub fn copy_vector(source: &Vec<f64>, destination: &mut Vec<f64>) {
    //for i in 0..source.len() {
    //    destination[i]=source[i];
    //}
    if source.len() != destination.len() {
        panic!("vectors have not the same length !!!");
    } else {
        destination[..source.len()].clone_from_slice(&source[..source.len()]);
    }
}

pub fn copy_solution(source: &Genome, destination: &mut Genome, dim: usize) {
    destination.id = source.id;
    destination.fitness = source.fitness;
    destination.genes[..dim].clone_from_slice(&source.genes[..dim]);
}

pub fn copy_vector2genome(source: &Vec<f64>, destination: &mut Genome) {
    for i in 0..source.len() {
        destination.genes[i] = source[i];
    }
}

///
/// Compute the euclidean distance between 2 solutions
///
#[allow(dead_code)]
pub fn euclidian_dist(p1: &Genome, p2: &Genome) -> f64 {
    let sum = p1
        .genes
        .iter()
        .zip(p2.genes.iter())
        .fold(0.0f64, |acc, (a, b)| acc + f64::powi(a - b, 2));
    f64::sqrt(sum)
}

/*
///
/// The S-Shape-V2() function is used to perform a binary optimization.
/// S-Shape-V2(x) = 1/(1+e^(-x))
///
#[cfg(feature = "binary")]
pub fn s_shape_v2(solution: &mut Genome, rng: &mut ThreadRng) {
    //V(i,:)=abs((2/pi)*atan((pi/2)*DeltaC(i,:)));
    let d: usize = solution.get_dimensions();

    let mut t_vec: Vec<f64> = vec![0.0; d];

    for j in 0..d {
        t_vec[j] = 1.0 / (1.0 + f64::exp(-1.0 * solution.genes[j]));
    }

    let intervall_01: Uniform<f64> = Uniform::from(0.0..=1.0);

    for j in 0..d {
        if intervall_01.sample(rng) < t_vec[j] {
            solution.genes[j] = 1.0;
        } else {
            solution.genes[j] = 0.0;
        }
    }
}
*/

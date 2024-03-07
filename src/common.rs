extern crate rand;
use rand::distributions::Uniform;
use rand::distributions::Distribution;

use crate::core::genome::Genome;


#[allow(dead_code)]
pub fn copy_matrix(source : &[Genome], destination : &mut Vec<Vec<f64>>) {
        
    let ni = source.len();
    let nj = source[0].get_dimensions();

    for i in 0..ni {
        for j in 0..nj {
            destination[i][j] =source[i].genes[j];
        }
     }
}

pub fn randomize(randvect : &mut Vec<f64>) {    
    let between = Uniform::from(0.0..=1.0);
    let mut rng = rand::thread_rng();
            
    for i in 0..randvect.len() {
        randvect[i]=between.sample(&mut rng);
    }
}

pub fn copy_vector(source : & Vec<f64>, destination : &mut Vec<f64>, dim : usize){
     
    //for i in 0..source.len() {
    //    destination[i]=source[i];
    //}
    destination[..dim].clone_from_slice(&source[..dim]);    


}

pub fn copy_vector2genome(source : & Vec<f64>, destination : &mut Genome){
    for i in 0..source.len() {
        destination.genes[i]=source[i];
    }
}


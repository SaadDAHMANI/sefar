extern crate rand;
use rand::distributions::Uniform;
use rand::distributions::Distribution;
use rayon::vec;

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


pub fn randi(imin: usize, imax : usize, cols : usize)->Vec<usize> {

    let mut vector : Vec<usize> = vec![0; cols];

    println!("i_min {}, imax : {}", imin, imax);
    if imin == imax {
        for j in 0..cols{
            vector[j] = imin;
        }    
    }
    else {
        let between = Uniform::from(imin..imax);
        let mut rng = rand::thread_rng();
                
        for j in 0..cols{
            vector[j] = between.sample(&mut rng);
        }      
    }
    
    vector
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

///
/// Compute the euclidean distance between 2 solutions
/// 
pub fn euclidian_dist(p1 : &Genome, p2 : &Genome)-> f64 {
   
    let sum  = p1.genes.iter().zip(p2.genes.iter()).fold(0.0f64, |acc, (a,b)| acc+f64::powi(a-b, 2) );
    
    f64::sqrt(sum)
        
    
}

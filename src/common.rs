extern crate rand;
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




extern crate rand;
use rand::distributions::Uniform;
use rand::distributions::Distribution;

use crate::core::parameters::Parameters;
use crate::core::genome::Genome;



pub fn initialize<P: Parameters>(params: &P)-> Vec<Genome>{

    let n: usize = params.get_population_size();
    let dim: usize = params.get_dimensions();
    let lb = params.get_lower_bounds();
    let ub = params.get_upper_bounds();
    
    let mut positions : Vec<Genome> = Vec::with_capacity(n);

    let intervall01 = Uniform::from(0.0f64..=1.0f64);
    let mut rng = rand::thread_rng();              
  
    for i in 0..n{
           let mut sln = Genome::new(i, dim); 
            
           for  j in 0..dim {   
            //  positions[i][j]= intervall01.sample(&mut rng)*(ub-lb)+lb;                         
            sln.genes[j]= intervall01.sample(&mut rng)*(ub[j]-lb[j])+lb[j];   
          }
          positions.push(sln);
    }        
    positions
}

#[allow(dead_code)]
pub fn copy_genome(source : &Genome, destination : &mut Genome){
    destination.id= source.id;
    destination.fitness = source.fitness;

    for i in 0..source.get_dimensions() {
        destination.genes[i]=source.genes[i];
    }
}

pub fn randomize(randvect : &mut Vec<f64>) {    
    let between = Uniform::from(0.0..=1.0);
    let mut rng = rand::thread_rng();
            
    for i in 0..randvect.len() {
        randvect[i]=between.sample(&mut rng);
    }
}

pub fn copy_vector(source : & Vec<f64>, destination : &mut Vec<f64>){
    for i in 0..source.len() {
        destination[i]=source[i];
    }
}

pub fn copy_vector2genome(source : & Vec<f64>, destination : &mut Genome){
    for i in 0..source.len() {
        destination.genes[i]=source[i];
    }
}

pub fn copy_matrix(source : & Vec<Genome>, destination : &mut Vec<Vec<f64>>) {
        
    let ni = source.len();
    let nj = source[0].get_dimensions();

    for i in 0..ni {
        for j in 0..nj {
            destination[i][j] =source[i].genes[j];
        }
     }
}

pub fn check_parameters<P: Parameters>(params : &P)-> Result<(), String> {
    
    let mut errors : usize = 0;
    let mut msg : String = String::new();
    
    if params.get_population_size() == 0 {
        msg = String::from("population_size must be greater than 0!; \n");
        errors +=1;
    }

    if params.get_dimensions() ==0 {
        msg = format!("{} Search space dimensions must be greater than 0!; \n", msg);
        errors +=1;
    }

    if params.get_max_iterations() ==0 {
        msg = format!("{} Iterations count (max_iterations) must be greater than 0!; \n", msg);
        errors +=1;
    }

    if params.get_lower_bounds().len() == 0 {
        msg = format!("{} Lower_bounds length must be greater than 0!; \n", msg);
        errors +=1;
    }

    if params.get_upper_bounds().len() == 0 {
        msg = format!("{} Upper_bounds length must be greater than 0!; \n", msg);
        errors +=1;
    }

    if params.get_lower_bounds().len() != params.get_upper_bounds().len() {
        msg = format!("{} Lower_bounds & Upper_bounds lengths must be equal!; \n", msg);
        errors +=1;
    }

    if params.get_lower_bounds().len() != params.get_dimensions() || params.get_upper_bounds().len() != params.get_dimensions() {
        msg = format!("{} Lower_bounds & Upper_bounds lengths must equal search space dimension!; \n", msg);
        errors +=1;
    }
   
    if errors > 0  {
        msg = format!("There are [{}] errors : \n {}", errors, msg.trim());       
        Err(msg)
    }
    else {
        Ok(())
    }      
}


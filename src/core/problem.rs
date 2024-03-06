//
// The Objectivefunction can be used with sequential algorithms (mutation of problem can be done).
//


pub trait Problem : Send + Sync + Clone { //+ 'static {
    fn objectivefunction(&mut self, genome : &[f64])->f64 {         
        genome.iter().fold(0.0f64, |sum, x| sum +x)
    }    
}

//pub trait ParaObjectivefunction : Send + Sync + 'static {
//    fn evaluate(&mut self, genome : &Vec<f64>)->f64 {         
//        genome.iter().fold(0.0f64, |sum, x| sum +x)
//    }    
//}

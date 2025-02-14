//use std::fmt::Display;
use crate::core::genome::Genome;
use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Duration;

#[derive(Debug, Clone)]
pub struct OptimizationResult {
    pub best_genome: Option<Genome>,
    pub best_fitness: Option<f64>,
    pub convergence_trend: Option<Vec<f64>>,
    pub computation_time: Option<Duration>,
    pub err_report: Option<String>,
}

#[allow(dead_code)]
impl OptimizationResult {
    pub fn get_none(msg: String) -> OptimizationResult {
        OptimizationResult {
            best_genome: None,
            best_fitness: None,
            convergence_trend: None,
            computation_time: None,
            err_report: Some(msg),
        }
    }
}

impl OptimizationResult {
    ///
    /// Save OptimizationResult instance in CSV file.
    ///
    #[allow(dead_code)]
    pub fn save(&self, header: Option<&str>, path: impl AsRef<Path>) -> std::io::Result<()> {
        let file = File::create(path.as_ref())?;

        let mut writer = BufWriter::new(file);
        match header {
            None => {}
            Some(txt) => {
                write!(writer, "{txt}\n")?;
            }
        }

        match &self.convergence_trend {
            None => {
                write!(writer, "Error-Report: {:?}\n", self.err_report)?;
            }
            Some(vector) => {
                write!(writer, "Best-fitness: {:?}\n", self.best_fitness)?;
                write!(writer, "Best-Solution: {:?}\n", self.best_genome)?;
                write!(writer, "Computation-Time: {:?}\n", self.computation_time)?;
                write!(writer, "Best-fitness history: \n")?;
                for x in vector.iter() {
                    write!(writer, "{x}\n")?;
                }
            }
        }

        writer.flush()?;

        Ok(())
    }
}

impl Display for OptimizationResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Best-fitness : {:?}; Best-solution : {:?}; Time : {:?}; Err-report: {:?}",
            self.best_fitness, self.best_genome, self.computation_time, self.err_report
        )
    }
}

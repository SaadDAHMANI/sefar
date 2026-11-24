//use std::fmt::Display;
use crate::core::genome::Genome;
use crate::core::OptError;

use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Duration;

#[derive(Debug, Clone)]
pub struct OptimizationResult {
    best_genome: Option<Genome>,
    best_fitness: Option<f64>,
    convergence_trend: Option<Vec<f64>>,
    computation_time: Option<Duration>,
    err_report: Option<OptError>,

    /// store the real number of parallel threads.
    #[cfg(feature = "parallel")]
    used_num_threads: usize,
}

#[allow(dead_code)]
impl OptimizationResult {
    pub fn get_empty(opt_error: Option<OptError>) -> OptimizationResult {
        OptimizationResult {
            best_genome: None,
            best_fitness: None,
            convergence_trend: None,
            computation_time: None,
            err_report: opt_error,

            #[cfg(feature = "parallel")]
            num_threads: 0,
        }
    }
}

impl OptimizationResult {
    #[cfg(not(feature = "parallel"))]
    pub fn new(
        best_genome: Option<Genome>,
        best_fitness: Option<f64>,
        convergence_trend: Option<Vec<f64>>,
        computation_time: Option<Duration>,
        err_report: Option<OptError>,
    ) -> Self {
        Self {
            best_genome,
            best_fitness,
            convergence_trend,
            computation_time,
            err_report,
        }
    }

    #[cfg(feature = "parallel")]
    pub fn new(
        best_genome: Option<Genome>,
        best_fitness: Option<f64>,
        convergence_trend: Option<Vec<f64>>,
        computation_time: Option<Duration>,
        err_report: Option<OptError>,
        num_threads: usize,
    ) -> Self {
        Self {
            best_genome,
            best_fitness,
            convergence_trend,
            computation_time,
            err_report,
            num_threads,
        }
    }

    pub fn best_genome(&self) -> Option<&Genome> {
        self.best_genome.as_ref()
    }

    pub fn best_fitness(&self) -> Option<f64> {
        self.best_fitness
    }

    pub fn convergence_trend(&self) -> Option<&Vec<f64>> {
        self.convergence_trend.as_ref()
    }

    pub fn computation_time(&self) -> Option<Duration> {
        self.computation_time
    }

    pub fn err_report(&self) -> Option<&OptError> {
        self.err_report.as_ref()
    }

    /// return the real number of parallel threads.
    #[cfg(feature = "parallel")]
    pub fn num_threads(&self) -> usize {
        self.num_threads
    }

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

                #[cfg(feature = "parallel")]
                write!(writer, "Number of parallel threads: {}\n", self.num_threads)?;

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

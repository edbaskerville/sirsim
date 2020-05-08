use std::path::{PathBuf, Path};
use std::fs::File;
use std::io::Read;

use sirtools::spec::*;

use serde::{Deserialize, Serialize};
use sirtools::stan::StanModel;

fn main() {
    let result = run();
    eprintln!("{:?}", result);
    
    // Write result to stdout
}

#[derive(Debug, Serialize, Deserialize)]
enum Error {
    InvalidInputPath(String),
    InvalidInputFile(String),
    InputReadFailure,
    InvalidJson(JsonError),
}

#[derive(Debug, Serialize, Deserialize)]
struct JsonError {
    description: String,
    category: String,
    line: usize,
    column: usize,
}

impl From<serde_json::error::Error> for Error {
    fn from(e: serde_json::error::Error) -> Self {
        Self::InvalidJson(JsonError {
            description: format!("{}", e),
            category: format!("{:?}", e.classify()),
            line: e.line(),
            column: e.column(),
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct InputData {
    structure: ModelStructure,
    config: ModelConfig,
}

#[derive(Debug, Serialize, Deserialize)]
struct OutputData {
    stan_code: String,
}

fn run() -> Result<OutputData, Error> {
    // Read config from file specified in first command-line argument or from stdin
    let args: Vec<String> = std::env::args().collect();
    let json_data = if args.len() > 1 {
        read_data_from_file(&args[1])?
    }
    else {
        read_data_from_stdin()?
    };
    let input_data: InputData = serde_json::from_str(&json_data)?;
    
    // Generate Stan code
    let stan_code = StanModel::new(
        input_data.structure, input_data.config
    ).generate_stan_code();
    eprintln!("{}", stan_code);
    
    Ok(OutputData {
        stan_code
    })
}

fn read_data_from_file(path_str: &str) -> Result<String, Error> {
    let path = Path::new(path_str).canonicalize().map_err(
        |_| Error::InvalidInputPath(path_str.into())
    )?;
    let mut file = File::open(&path).map_err(
        |_| Error::InvalidInputFile(path_str.into())
    )?;
    let mut data = String::new();
    file.read_to_string(&mut data).map_err(
        |_| Error::InputReadFailure
    )?;
    Ok(data)
}

fn read_data_from_stdin() -> Result<String, Error> {
    let mut data = String::new();
    std::io::stdin().read_to_string(&mut data).map_err(
        |_| Error::InputReadFailure
    )?;
    Ok(data)
}


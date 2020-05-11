use serde::{Serialize, Deserialize};

#[derive(Debug, Serialize, Deserialize)]
pub enum Error {
    InvalidInputPath(String),
    InvalidInputFile(String),
    InputReadFailure,
    InvalidJson(JsonError),
}

#[derive(Debug, Serialize, Deserialize)]
pub struct JsonError {
    pub description: String,
    pub category: String,
    pub line: usize,
    pub column: usize,
}

use std::path::{PathBuf, Path};
use std::fs::File;
use std::io::Read;

use sirtools::spec::*;

fn main() {
    let result = sirtools::stan::run();
    
    let output_json = serde_json::to_string_pretty(&result).unwrap();
    
    println!("{}", output_json);
}

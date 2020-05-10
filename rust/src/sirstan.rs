use std::path::{PathBuf, Path};
use std::fs::File;
use std::io::Read;

use sirtools::spec::*;

fn main() {
    let result = sirtools::stan::run();
    eprintln!("{:?}", result);
    
    // Write result to stdout
}

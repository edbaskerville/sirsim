use crate::errors::*;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;
use rusqlite::types::ValueRef;
use std::iter::FromIterator;

pub fn write_data_to_file(path_str: &str, data: &str) {
    let mut file = File::create(path_str).unwrap();
    file.write_all(data.as_bytes()).unwrap();
}

pub fn read_data_from_stdin() -> Result<String, Error> {
    let mut data = String::new();
    std::io::stdin().read_to_string(&mut data).map_err(
        |_| Error::InputReadFailure
    )?;
    Ok(data)
}

pub fn read_data_from_file(path_str: &str) -> Result<String, Error> {
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

pub fn db_table_to_json_object(
    conn: &rusqlite::Connection, table_name: &str,
    column_names: &Vec<&str>
) -> serde_json::Value {
    let col_values_pairs: Vec<(String, serde_json::Value)> = column_names.iter().map(|c| {
        let values: Vec<serde_json::Value> = conn.prepare(
            &format!("SELECT {} FROM {};", c, table_name)
        ).unwrap().query_map(rusqlite::params![], |row| {
            Ok(
                match row.get_raw(0) {
                    ValueRef::Null => {
                        serde_json::Value::Null
                    },
                    ValueRef::Integer(val) => {
                        val.into()
                    },
                    ValueRef::Real(val) => {
                        val.into()
                    },
                    ValueRef::Text(val) => {
                        serde_json::Value::String(String::from_utf8(val.to_vec()).unwrap())
                    },
                    ValueRef::Blob(_) => {
                        unimplemented!()
                    },
                }
            )
        }).unwrap().map(|r| r.unwrap()).collect::<Vec<_>>();
        
        (String::from(*c), serde_json::Value::Array(values))
    }).collect();
    
    let map = serde_json::Map::from_iter(col_values_pairs);
    
    serde_json::Value::Object(map)
}

#![allow(non_snake_case)]

use sirtools::ibm::*;
use std::time::Instant;

use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Read;
use std::path::{PathBuf, Path};
use std::collections::HashMap;
use std::f64::INFINITY;

use sirtools::util::*;
use sirtools::errors::*;
use std::iter::FromIterator;

#[derive(Serialize, Deserialize)]
struct Config {
    rng_seed: Option<u32>,
    output_path: Option<String>,
    write_to_stdout: Option<bool>,
    
    record_all_events: bool,
    
    t_initial: Option<f64>,
    t_final: Option<f64>,
    
    n_ageclasses: usize,
    
    susceptible_state: String,
    initial_infected_state: String,
    final_states: Vec<String>,
    
    infected_states: Vec<StateConfig>,
    
    contact_parameters: Vec<ContactParameters>,
    
    initial_counts: HashMap<String, Vec<usize>>,
}

#[derive(Serialize, Deserialize)]
struct ContactParameters {
    beta: f64,
    C: Vec<Vec<f64>>,
    t_end: Option<f64>,
}

#[derive(Serialize, Deserialize)]
struct StateConfig {
    name: String,
    infectious: bool,
    mean_duration: f64,
    gamma_shape: f64,
    next_states: Vec<String>,
    probabilities: Option<Vec<Vec<f64>>>,
}

fn main() -> Result<(), Error> {
    // Read JSON data from file specified in first command-line argument or from stdin
    let args: Vec<String> = std::env::args().collect();
    let json_data = if args.len() > 1 {
        read_data_from_file(&args[1])?
    }
    else {
        read_data_from_stdin()?
    };
    
    // Read config from JSON data
    let config: Config = serde_json::from_str(&json_data).unwrap();
    let (
        states,
        susceptible_state_id,
        initial_infected_state_id,
        initial_counts
    ) = parse_states(&config);
    let (t_change, beta_t, C_t) = parse_contact_parameters(config.contact_parameters);
    
    // If we were given a config file, use its parent as our working directory
    if args.len() > 1 {
        std::env::set_current_dir(&Path::new(&args[1]).parent().unwrap()).unwrap();
    }
    
    // Write to DB file specified in config file
    // (or use in-memory database if not specified)
    let mut db_connection = match config.output_path {
        Some(output_path) => {
            let db_path: PathBuf = output_path.into();
            assert!(!db_path.exists());
            rusqlite::Connection::open(db_path).unwrap()
        },
        None => {
            assert!(!config.record_all_events);
            rusqlite::Connection::open_in_memory().unwrap()
        }
    };
    
    let mut sim = {
        let mut db_transaction = db_connection.transaction().unwrap();
        let sim = Simulation::new(
            config.n_ageclasses,
            states,
            susceptible_state_id,
            initial_infected_state_id,
            t_change,
            beta_t,
            C_t,
            initial_counts,
            config.rng_seed,
            &mut db_transaction,
            config.record_all_events,
        );
        sim.write_counts(&mut db_transaction);
        db_transaction.commit().unwrap();
        sim
    };
    
    let start = Instant::now();
    let t_final = config.t_final.unwrap_or(INFINITY);
    println!("t = {}", sim.t);
    let mut done = false;
    while sim.t < t_final && !done {
        let mut db_transaction = db_connection.transaction().unwrap();
        done = sim.simulate(sim.t + 1.0, &mut db_transaction, config.record_all_events);
        sim.write_counts(&mut db_transaction);
        db_transaction.commit().unwrap();
        println!("t = {}", sim.t);
    }
    eprintln!("elapsed time: {} s", start.elapsed().as_secs_f64());
    
    eprintln!("...done.");
    
    if config.write_to_stdout.unwrap_or(false) {
        eprintln!("Writing DB to stdout in JSON format...");
        
        let db_json_data = serde_json::Map::from_iter(vec![
            ("Meta", vec!["key", "value"]),
            ("Counts", vec!["time", "state", "ageclass", "count"]),
            ("RtSufficientStatistics", vec!["time_discrete", "n_primary", "n_secondary"]),
        ].iter().map(|(table_name, col_names)| {
            (
                String::from(*table_name),
                db_table_to_json_object(&db_connection, table_name, col_names)
            )
        }).collect::<Vec<_>>());
        
        println!("{}", serde_json::to_string_pretty(&db_json_data).unwrap());
    }
    
    Ok(())
}

fn parse_states(config: &Config) -> (Vec<State>, usize, usize, Counts) {
    let mut states = Vec::new();
    let mut n_states: usize = 0;
    let mut name_id_map = HashMap::new();
    
    // Susceptible state
    let susceptible_state_id = n_states;
    name_id_map.insert(config.susceptible_state.clone(), susceptible_state_id);
    states.push(State::new_susceptible(susceptible_state_id, config.susceptible_state.clone()));
    n_states += 1;
    
    // Final states
    for name in &config.final_states {
        name_id_map.insert(name.clone(), n_states);
        states.push(State::new_final(n_states, name.clone()));
        n_states += 1;
    }
    
    // Infected states, without transitions
    for state_config in &config.infected_states {
        name_id_map.insert(state_config.name.clone(), n_states);
        states.push(State::new_infected(
            n_states, state_config.name.clone(),
        ));
        n_states += 1;
    }
    
    // Add transitions, resolving to state IDs
    for state_config in &config.infected_states {
        let id = name_id_map[&state_config.name];
        
        let next_state_ids = state_config.next_states.iter().map(
            |name| name_id_map[name]
        ).collect();
        
        let transition_probabilities = match &state_config.probabilities {
            Some(probabilities) => {
                probabilities.clone()
            },
            None => {
                std::iter::repeat(Vec::new()).take(config.n_ageclasses).collect()
            }
        };
        
        let transition_cdfs = transition_probabilities.iter().map(|tprobs| {
            cumulative_sum(tprobs)
        }).collect();
        
        states[id].detail = StateDetail::Infected(Some(InfectedState {
            infectious: state_config.infectious,
            mean_duration: state_config.mean_duration,
            gamma_shape: state_config.gamma_shape,
            next_state_ids: next_state_ids,
            transition_cdfs: transition_cdfs,
        }))
    }
    
    let initial_infected_state_id = name_id_map[&config.initial_infected_state];
    
    let initial_counts = parse_initial_counts(
        states.len(), config.n_ageclasses,
        &name_id_map, &config.initial_counts
    );
    
    (states, susceptible_state_id, initial_infected_state_id, initial_counts)
}

fn parse_initial_counts(
    n_states: usize, n_ageclasses: usize,
    name_id_map: &HashMap<String, usize>, counts_raw: &HashMap<String, Vec<usize>>
) -> Counts {
    let mut counts = Counts::new(n_states, n_ageclasses);
    for (state_name, counts_for_state) in counts_raw {
        for i in 0..n_ageclasses {
            counts.increment(name_id_map[state_name], i, counts_for_state[i]);
        }
    }
    counts
}

fn parse_contact_parameters(cp_vec: Vec<ContactParameters>) -> (Vec<f64>, Vec<f64>, Vec<Vec<Vec<f64>>>) {
    let mut t_change = Vec::new();
    let mut beta_t = Vec::new();
    let mut C_t = Vec::new();
    
    for i in 0..cp_vec.len() {
        beta_t.push(cp_vec[i].beta);
        C_t.push(cp_vec[i].C.clone());
        
        if let Some(t_end) = cp_vec[i].t_end {
            assert!(i < cp_vec.len() - 1);
            t_change.push(t_end);
        }
        else {
            assert!(i == cp_vec.len() - 1);
        }
    }
    
    (t_change, beta_t, C_t)
}

use serde::{Serialize, Deserialize};
use std::collections::HashMap;

/// The structure of a compartmental epidemiological model,
/// including infection states and auxiliary variables used
/// to model observation delays.
#[derive(Debug, Serialize, Deserialize)]
pub struct ModelStructure {
    pub susceptible_state: String,
    pub states: Vec<State>,
    pub observation_variables: Vec<ObservationVariable>
}

/// A single state in the compartmental model.
#[derive(Debug, Serialize, Deserialize)]
pub struct State {
    pub name: String,
    pub infectious: bool,
    pub next_states: Option<Vec<String>>,
}

/// An auxiliary variable that accumulates delayed observations of
/// transitions between two states in the underlying infection process.
#[derive(Debug, Serialize, Deserialize)]
pub struct ObservationVariable {
    pub name: String,
    pub start_state: String,
    pub end_state: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ModelConfig {
    pub infer_process_delays: HashMap<String, bool>,
    pub infer_observation_delays: HashMap<String, bool>,
    pub observation_distributions: HashMap<String, ObservationDistribution>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(tag = "distribution", content = "infer")]
pub enum ObservationDistribution {
    Normal(NormalObservationParameters),
    Poisson(PoissonObservationParameters),
    NegativeBinomial(NegativeBinomialObservationParameters),
    Binomial(BinomialObservationParameters),
    BetaBinomial(BetaBinomialObservationParameters)
}

#[derive(Debug, Serialize, Deserialize)]
pub struct NormalObservationParameters {
    pub mean_fraction: bool,
    pub standard_deviation: bool
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PoissonObservationParameters {
    pub mean_fraction: bool,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct NegativeBinomialObservationParameters {
    pub mean_fraction: bool,
    pub dispersion: bool,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BinomialObservationParameters {
    pub probability: bool,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BetaBinomialObservationParameters {
    pub probability: bool,
    pub dispersion: bool,
}

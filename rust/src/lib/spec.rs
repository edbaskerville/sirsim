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
    pub population_size: f64,
    pub transmission_covariates: Vec<Covariate>,
    pub process_delays: HashMap<String, DelayConfig>,
    pub observations: HashMap<String, ObservationConfig>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Covariate {
    pub name: String,
    pub values: Vec<(f64, f64)>, // (time, value)
}

#[derive(Debug, Serialize, Deserialize)]
pub struct DelayConfig {
    pub mean_duration: Parameter,
    pub gamma_shape: Parameter,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ObservationConfig {
    pub delay: DelayConfig,
    pub distribution: ObservationDistribution,
    pub values: Vec<(f64, f64, f64)>, // (start_time, end_time, value)
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(tag = "tag", content = "content")]
pub enum Parameter {
    Fixed(f64),
    Inferred(ParameterDistribution),
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(tag = "tag", content = "content")]
pub enum ParameterDistribution {
    RawCode(String),
    Structured(StructuredParameterDistribution),
}

#[derive(Debug, Serialize, Deserialize)]
pub struct StructuredParameterDistribution {
    pub type_: String,
    pub parameters: Vec<Parameter>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(tag = "tag", content = "content")]
pub enum ObservationDistribution {
    Normal(NormalObservationParameters),
    Poisson(PoissonObservationParameters),
    NegativeBinomial(NegativeBinomialObservationParameters),
    Binomial(BinomialObservationParameters),
    BetaBinomial(BetaBinomialObservationParameters)
}

#[derive(Debug, Serialize, Deserialize)]
pub struct NormalObservationParameters {
    pub mean_fraction: Parameter,
    pub standard_deviation: Parameter
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PoissonObservationParameters {
    pub mean_fraction: Parameter,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct NegativeBinomialObservationParameters {
    pub mean_fraction: Parameter,
    pub dispersion: Parameter,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BinomialObservationParameters {
    pub probability: Parameter,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BetaBinomialObservationParameters {
    pub probability: Parameter,
    pub dispersion: Parameter,
}

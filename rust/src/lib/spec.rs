use serde::{Serialize, Deserialize};
use std::collections::HashMap;

/// The structure of a compartmental epidemiological model,
/// including infection states and auxiliary variables used
/// to model observation delays.
#[derive(Debug, Serialize, Deserialize)]
pub struct ModelStructure {
    susceptible_state: String,
    states: Vec<State>,
    observation_variables: Vec<ObservationVariable>
}

/// A single state in the compartmental model.
#[derive(Debug, Serialize, Deserialize)]
pub struct State {
    name: String,
    infectious: bool,
    next_states: Option<Vec<String>>,
}

/// An auxiliary variable that accumulates delayed observations of
/// transitions between two states in the underlying infection process.
#[derive(Debug, Serialize, Deserialize)]
pub struct ObservationVariable {
    name: String,
    start_state: String,
    end_state: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ModelConfig {
    population_size: f64,
    transmission_covariates: Vec<Covariate>,
    process_delays: HashMap<String, DelayConfig>,
    observations: HashMap<String, ObservationConfig>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Covariate {
    name: String,
    values: Vec<(f64, f64)>, // (time, value)
}

#[derive(Debug, Serialize, Deserialize)]
pub struct DelayConfig {
    mean_duration: Parameter,
    gamma_shape: Parameter,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ObservationConfig {
    delay: DelayConfig,
    distribution: ObservationDistribution,
    values: Vec<(f64, f64, f64)>, // (start_time, end_time, value)
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
    type_: String,
    parameters: Vec<Parameter>,
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
    mean_fraction: Parameter,
    standard_deviation: Parameter
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PoissonObservationParameters {
    mean_fraction: Parameter,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct NegativeBinomialObservationParameters {
    mean_fraction: Parameter,
    dispersion: Parameter,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BinomialObservationParameters {
    probability: Parameter,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BetaBinomialObservationParameters {
    probability: Parameter,
    dispersion: Parameter,
}

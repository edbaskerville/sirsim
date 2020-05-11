use crate::spec::*;
use std::collections::HashMap;
use std::convert::{TryFrom, TryInto};


use serde::{Deserialize, Serialize};
use std::path::Path;
use std::fs::File;
use std::io::{Read, Write};

use indoc::indoc;
use unindent::unindent;

pub struct StanModel {
    structure: ModelStructure,
    config: ModelConfig,
}

impl StanModel {
    pub fn new(structure: ModelStructure, config: ModelConfig) -> Self {
        Self { structure, config }
    }
    
    pub fn generate_stan_code(&self) -> String {
        generate_stan_code(
            self.functions_code(),
            self.data_code(),
            self.transformed_data_code(),
            self.parameters_code(),
            self.transformed_parameters_code(),
            self.model_code(),
            self.generated_quantities_code(),
        )
    }
    
    pub fn functions_code(&self) -> Vec<String> {
        vec![
            self.ddt_function_code()
        ]
    }
    
    fn n_states(&self) -> usize {
        self.structure.states.len()
    }
    
    fn susceptible_state(&self) -> String {
        self.structure.susceptible_state.clone()
    }
    
    fn ddt_function_code(&self) -> String {
        let declaration = indoc!("
            real[] ode_ddt(
              real t, real[] state,
              real[] params,
              real[] x_r, int[] x_i
            )"
        );
        
        let body_sections = vec![
            "int n_substates;".into(),
            self.gamma_shape_declarations(),
            self.standard_declarations(true, true, false),
            self.mean_duration_declarations(true, true, false),
            self.probability_declarations(false),
            self.extract_or_assign("x_i", self.x_i_entries(), true),
            self.extract_or_assign("x_r", self.x_r_entries(), true),
            self.extract_or_assign("params", self.params_entries(), true),
            self.ddt_computation(),
        ];
        
        format!(
            "{declaration} {{\n{body}\n}}",
            declaration = declaration,
            body = indent(&body_sections.join("\n\n"), 2, 1)
        )
    }
    
    fn gamma_distributed_states(&self) -> Vec<String> {
        let mut vars = Vec::new();
        for state_config in &self.structure.states {
            if state_config.name != self.susceptible_state() {
                if let Some(next_states) = &state_config.next_states {
                    if next_states.len() > 0 {
                        vars.push(state_config.name.clone());
                    }
                }
            }
        }
        vars
    }
    
    fn observation_variables(&self) -> Vec<String> {
        self.structure.observation_variables.iter().map(
            |x| x.name.clone()
        ).collect()
    }
    
    // Variables with delays
    fn variables_with_delays(&self, data: bool, params: bool) -> Vec<String> {
        let mut v = Vec::new();
        
        for var in self.gamma_distributed_states() {
            if self.config.infer_process_delays[&var] {
                if params {
                    v.push(var);
                }
            }
            else {
                if data {
                    v.push(var);
                }
            }
        }
    
        for var in self.observation_variables() {
            if self.config.infer_observation_delays[&var] {
                if params {
                    v.push(var);
                }
            }
            else {
                if data {
                    v.push(var);
                }
            }
        }
        
        v
    }
    
    // Gamma shape for variables with delays
    fn gamma_shape_declarations(&self) -> String {
        self.variables_with_delays(true, true).iter().map(|var| {
            format!("int {}_gamma_shape;", var)
        }).collect::<Vec<_>>().join("\n")
    }
    
    // Fixed parameter declarations
    fn standard_declarations(&self, data: bool, params: bool, include_bounds: bool) -> String {
        let mut v: Vec<String> = Vec::new();
        
        if data {
            v.push("N".into());
        }
        
        if params {
            v.push("b".into());
        }
        
        v.iter().map(|var| {
            format!(
                "real{} {};",
                if include_bounds { "<lower=0>" } else { "" },
                var
            )
        }).collect::<Vec<_>>().join("\n")
    }
    
    // Durations for variables with delays
    fn mean_duration_declarations(&self, data: bool, params: bool, include_bounds: bool) -> String {
        self.variables_with_delays(data, params).iter().map(|var| {
            format!(
                "real{} {}_mean_duration;",
                if include_bounds { "<lower=0>" } else { "" },
                var
            )
        }).collect::<Vec<_>>().join("\n")
    }
    
    fn transitions(&self) -> Vec<(String, String)> {
        let mut transitions = Vec::new();
        for state_config in &self.structure.states {
            if let Some(next_states) = &state_config.next_states {
                for next_state in next_states {
                    transitions.push((state_config.name.clone(), next_state.clone()));
                }
            }
        }
        transitions
    }
    
    fn transitions_with_probabilities(&self) -> Vec<(String, String)> {
        let mut transitions = Vec::new();
        for state_config in &self.structure.states {
            if let Some(next_states) = &state_config.next_states {
                if next_states.len() > 1 {
                    for next_state in next_states {
                        transitions.push((state_config.name.clone(), next_state.clone()));
                    }
                }
            }
        }
        transitions
    }
    
    fn probability_declarations(&self, include_bounds: bool) -> String {
        self.transitions_with_probabilities().iter().map(|(from, to)| {
            format!(
                "real{} p_{}_{};",
                if include_bounds { "<lower=0, upper=1>" } else { "" },
                from, to
            )
        }).collect::<Vec<_>>().join("\n")
    }
    
    fn extract_or_assign(&self, vec_name: &str, vars: Vec<String>, extract: bool) -> String {
        let vec_access = format!("{}[index]", vec_name);
        let increment: String = "index += 1;".into();
        
        let mut sections = Vec::new();
        
        sections.push("int index = 1;".into());
        
        for var in vars {
            sections.push(vec![
                format!(
                    "{} = {};",
                    if extract { &var } else { &vec_access },
                    if extract { &vec_access } else { &var },
                ),
                increment.clone(),
            ].join("\n"));
        }
        
        format_block(sections.join("\n\n"))
    }
    
    fn x_i_entries(&self) -> Vec<String> {
        let mut v = Vec::new();
        
        v.push("n_substates".into());
        
        for var in self.variables_with_delays(true, true) {
            v.push(format!("{}_gamma_shape", var));
        }
        
        v
    }
    
    fn x_r_entries(&self) -> Vec<String> {
        let mut v = Vec::new();
        
        v.push("N".into());
        
        for var in self.fixed_delay_vars() {
            v.push(format!("{}_mean_duration", var));
        }
    
        for (from, to) in self.transitions_with_probabilities() {
            v.push(format!("p_{}_{}", from, to));
        }
        
        v
    }
    
    fn inferred_delay_vars(&self) -> Vec<String> {
        let mut v = Vec::new();
        
        for var in self.gamma_distributed_states() {
            if self.config.infer_process_delays[&var] {
                v.push(var);
            }
        }
        for var in self.observation_variables() {
            if self.config.infer_observation_delays[&var] {
                v.push(var);
            }
        }
        
        v
    }
    
    fn fixed_delay_vars(&self) -> Vec<String> {
        let mut v = Vec::new();
        
        for var in self.gamma_distributed_states() {
            if !self.config.infer_process_delays[&var] {
                v.push(var);
            }
        }
        for var in self.observation_variables() {
            if !self.config.infer_observation_delays[&var] {
                v.push(var);
            }
        }
        
        v
    }
    
    fn params_entries(&self) -> Vec<String> {
        let mut v = Vec::new();
        
        v.push("b".into());
        for var in self.inferred_delay_vars() {
                v.push(format!("{}_mean_duration", var));
        }
        
        v
    }
    
    fn ddt_computation(&self) -> String {
        let sections = vec![
            self.ddt_computation_declarations(),
            self.ddt_extract_state(),
            self.ddt_transmission(),
            self.ddt_transitions_between_states(),
            self.ddt_transitions_within_states(),
            self.ddt_state_changes(),
            self.ddt_obs_var_changes(),
            self.ddt_assembly(),
            "return ddt;".into(),
        ];
        
        format_block(sections.join("\n\n"))
    }
    
    fn ddt_computation_declarations(&self) -> String {
        let sections = vec![
            "real ddt[n_substates];".into(),
            self.ddt_state_declarations(),
            self.ddt_transitions_between_states_declarations(),
            self.ddt_transitions_within_states_declarations(),
            self.ddt_state_changes_declarations(),
        ];
        
        sections.join("\n\n")
    }
    
    fn state_isgamma_pairs(&self) -> Vec<(String, bool)> {
        self.structure.states.iter().filter_map(|state| {
            if state.name.eq(&self.susceptible_state()) {
                None
            }
            else {
                Some(
                    if let Some(next_states) = &state.next_states {
                        (state.name.clone(), true)
                    }
                    else {
                        (state.name.clone(), false)
                    }
                )
            }
        }).collect()
    }
    
    fn ddt_state_declarations(&self) -> String {
        let mut lines = Vec::new();
        
        // States involved in epidemic process
        lines.push(format!("real {};", self.susceptible_state()));
        for (name, has_delay) in &self.state_isgamma_pairs() {
            if *has_delay {
                lines.push(
                    format!("real {0}[{0}_gamma_shape];", name)
                );
            }
            else {
                // Final states
                lines.push(
                    format!("real {};", name)
                );
            }
        }
        
        // Observation variables
        for name in &self.observation_variables() {
            lines.push(
                format!("real {0}[{0}_gamma_shape + 1];", name)
            );
        }
        
        lines.join("\n")
    }
    
    fn ddt_transitions_between_states_declarations(&self) -> String {
        self.transitions().iter().map(|(from, to)| {
            format!("real d_{}_{};", from, to)
        }).collect::<Vec<_>>().join("\n")
    }
    
    fn ddt_transitions_within_states_declarations(&self) -> String {
        let mut lines = Vec::new();
        
        for name in self.gamma_distributed_states() {
            lines.push(
                format!("real d_{0}_{0}[{0}_gamma_shape - 1];", name)
            );
        }
        
        for name in self.observation_variables() {
            lines.push(
                format!("real d_{0}_{0}[{0}_gamma_shape];", name)
            );
        }
        
        lines.join("\n")
    }
    
    fn ddt_state_changes_declarations(&self) -> String {
        let mut lines = Vec::new();
        
        for (name, has_delay) in self.state_isgamma_pairs() {
            if has_delay {
                lines.push(
                    format!("real d_{0}[{0}_gamma_shape];", name)
                );
            }
            else {
                lines.push(
                    format!("real d_{};", name)
                );
            }
        }
        
        for name in self.observation_variables() {
            lines.push(
                format!("real d_{0}[{0}_gamma_shape + 1];", name)
            );
        }
        
        lines.join("\n")
    }
    
    fn ddt_extract_state(&self) -> String {
        let mut sections = Vec::new();
        
        sections.push("int index = 1;".into());
        
        for (name, isgamma) in self.state_isgamma_pairs() {
            sections.push({
                if isgamma {
                    vec![
                        format!("{0} = state[index:(index + {0}_gamma_shape - 1)];", name),
                        format!("index += {}_gamma_shape;", name),
                    ]
                }
                else {
                    vec![
                        format!("{0} = state[index];", name),
                        "index += 1;".into(),
                    ]
                }
            }.join("\n"))
        }
    
        for name in self.observation_variables() {
            sections.push(vec![
                format!("{0} = state[index:(index + {0}_gamma_shape)];", name),
                format!("index += {}_gamma_shape + 1;", name),
            ].join("\n"));
        }
        
        sections.push({
            let mut lines = Vec::new();
            
            lines.push(format!("{} = N", self.susceptible_state()));
            
            for (name, isgamma) in self.state_isgamma_pairs() {
                if isgamma {
                    lines.push(format!("  - sum({})", name));
                }
                else {
                    lines.push(format!("  - {}", name));
                }
            }
            
            lines.push(";".into());
            
            lines.join("\n")
        });
        
        format_block(sections.join("\n\n"))
    }
    
    fn ddt_transmission(&self) -> String {
        let S = self.susceptible_state();
        let next_state = self.next_states(&S)[0].clone();
        
        let sum_infectious: Vec<String> = self.infectious_states().iter().map(|name| {
            format!("sum({})", name)
        }).collect();
        
        format_block(
            format!(
                "d_{S}_{next_state} = b * ({sum_infectious}) * {S} / N;",
                S = S, next_state = next_state, sum_infectious = sum_infectious.join(" + ")
            )
        )
    }
    
    fn next_states(&self, name: &String) -> Vec<String> {
        for state in &self.structure.states {
            if state.name.eq(name) {
                if let Some(next_states) = &state.next_states {
                    return next_states.clone()
                }
                else {
                    return vec![]
                }
            }
        }
        vec![]
    }
    
    fn infectious_states(&self) -> Vec<String> {
        self.structure.states.iter().filter_map(|state| {
            if state.infectious {
                Some(state.name.clone())
            }
            else {
                None
            }
        }).collect()
    }
    
    fn ddt_transitions_between_states(&self) -> String {
        let mut sections = Vec::new();
        
        for state in &self.structure.states {
            if state.name != self.susceptible_state() {
                if let Some(next_states) = &state.next_states {
                    if next_states.len() == 1 {
                        sections.push(format_block(format!(
                            "d_{0}_{1} = {0}[{0}_gamma_shape] * {0}_gamma_shape / {0}_mean_duration;",
                            state.name, next_states[0]
                        )));
                    }
                    else if next_states.len() > 1 {
                        let mut lines = Vec::new();
                        
                        lines.push(format!(
                            "real d_{0}_ = {0}[{0}_gamma_shape] * {0}_gamma_shape / {0}_mean_duration;",
                            state.name
                        ));
                        for next_state in next_states {
                            lines.push(format!(
                                "d_{0}_{1} = p_{0}_{1} * d_{0}_;",
                                state.name, next_state
                            ));
                        }
                        sections.push(format_block(lines.join("\n")));
                    }
                }
            }
        }
        
        sections.join("\n\n")
    }
    
    fn ddt_transitions_within_states(&self) -> String {
        let mut sections = Vec::new();
        
        for name in self.gamma_distributed_states() {
            sections.push(vec![
                format!("for(i in 1:({}_gamma_shape - 1)) {{", name),
                format!("  d_{0}_{0}[i] = {0}[i] * {0}_gamma_shape / {0}_mean_duration;", name),
                "}".into()
            ].join("\n"));
        }
        
        for name in self.observation_variables() {
            sections.push(vec![
                format!("for(i in 1:{}_gamma_shape) {{", name),
                format!("  d_{0}_{0}[i] = {0}[i] * {0}_gamma_shape / {0}_mean_duration;", name),
                "}".into()
            ].join("\n"));
        }
        
        sections.join("\n\n")
    }
    
    fn ddt_state_changes(&self) -> String {
        let mut sections = Vec::new();
        
        for (state, isgamma) in self.state_isgamma_pairs() {
            let prev_state = self.previous_state(&state);
            let next_states = self.next_states(&state);
            
            if isgamma {
                let subtract_next_states = next_states.iter().map(|next_state|
                    format!(" - d_{}_{}", state, next_state)
                ).collect::<Vec<_>>().join("");
                
                let if_gamma1 = format!(
                    "d_{0}[1] = d_{1}_{0}{2};",
                    state, prev_state, subtract_next_states
                );
                
                let else_gamma1 = vec![
                    format!("d_{0}[1] = d_{1}_{0} - d_{0}_{0}[1];", state, prev_state),
                    format!("for(i in 2:({}_gamma_shape - 1)) {{", state),
                        format!("  d_{0}[i] = d_{0}_{0}[i - 1] - d_{0}_{0}[i];", state),
                    "}".into(),
                    format!(
                        "d_{0}[{0}_gamma_shape] = d_{0}_{0}[{0}_gamma_shape - 1]{1};",
                        state, subtract_next_states
                    ),
                ].join("\n");
                
                sections.push(vec![
                    format!("if({}_gamma_shape == 1) {{", state),
                        indent(&if_gamma1, 2, 1),
                    "}".into(),
                    "else {".into(),
                        indent(&else_gamma1, 2, 1),
                    "}".into()
                ].join("\n"));
            }
            else {
                sections.push(format_block(
                    format!("d_{0} = d_{1}_{0};", state, prev_state)
                ));
            }
        }
        
        sections.join("\n\n")
    }
    
    fn previous_state(&self, name: &String) -> String {
        for state in &self.structure.states {
            if let Some(next_states) = &state.next_states {
                if next_states.contains(name) {
                    return state.name.clone();
                }
            }
        }
        panic!()
    }
    
    fn ddt_obs_var_changes(&self) -> String {
        let mut sections = Vec::new();
    
        for obs_var in &self.structure.observation_variables {
            let name = obs_var.name.clone();
            let start_state = obs_var.start_state.clone();
            let end_state = obs_var.end_state.clone();
            
            let if_gamma0 = format!(
                "d_{0}[1] = d_{1}_{2};",
                name, start_state, end_state
            );
        
            let else_gamma0 = vec![
                format!(
                    "d_{0}[1] = d_{1}_{2} - d_{0}_{0}[1];",
                    name, start_state, end_state
                ),
                format!("for(i in 2:{}_gamma_shape) {{", name),
                    format!("  d_{0}[i] = d_{0}_{0}[i - 1] - d_{0}_{0}[i];", name),
                "}".into(),
                format!(
                    "d_{0}[{0}_gamma_shape + 1] = d_{0}_{0}[{0}_gamma_shape];",
                    name
                ),
            ].join("\n");
        
            sections.push(vec![
                format!("if({}_gamma_shape == 0) {{", name),
                    indent(&if_gamma0, 2, 1),
                "}".into(),
                "else {".into(),
                    indent(&else_gamma0, 2, 1),
                "}".into()
            ].join("\n"));
        }
    
        sections.join("\n\n")
    }
    
    fn ddt_assembly(&self) -> String {
        let mut sections = Vec::new();
        
        sections.push("int index = 1;".into());
        
        for (state, isgamma) in self.state_isgamma_pairs() {
            if isgamma {
                sections.push(vec![
                    format!("ddt[index:(index + {0}_gamma_shape - 1)] = d_{0};", state),
                    format!("index += {}_gamma_shape;", state),
                ].join("\n"));
            }
            else {
                sections.push(vec![
                    format!("ddt[index] = d_{0};", state),
                    "index += 1;".into(),
                ].join("\n"));
            }
        }
        
        for obs_var in self.observation_variables() {
            sections.push(vec![
                format!("ddt[index:(index + {0}_gamma_shape)] = d_{0};", obs_var),
                format!("index += {}_gamma_shape + 1;", obs_var),
            ].join("\n"));
        }
        
        format_block(sections.join("\n\n"))
    }
    
    pub fn data_code(&self) -> String {
        vec![
            self.gamma_shape_declarations(),
            self.standard_declarations(true, false, true),
            self.mean_duration_declarations(true, false, true),
            self.probability_declarations(true),
            self.ic_declarations(),
            self.time_declarations(),
        ].join("\n\n")
    }
    
    fn ic_declarations(&self) -> String {
        let mut lines = Vec::new();
        
        for state in &self.structure.states {
            if state.name != self.susceptible_state() {
                lines.push(format!("real<lower=0> {}_init;", state.name));
            }
        }
        
        lines.join("\n")
    }
    
    fn time_declarations(&self) -> String {
        vec![
            String::from("int n_times;"),
            "real start_time;".into(),
            "real times[n_times];".into(),
        ].join("\n")
    }
    
    fn transformed_data_code(&self) -> String {
        vec![
            self.transformed_data_declarations(),
            self.extract_or_assign("x_i", self.x_i_entries(), false),
            self.extract_or_assign("x_r", self.x_r_entries(), false),
        ].join("\n\n")
    }
    
    fn transformed_data_declarations(&self) -> String {
        vec![
            format!(
                "int n_delay_vars = {};",
                self.gamma_distributed_states().len() + self.observation_variables().len()
            ),
            format!(
                "int n_inferred_delays = {};", self.inferred_delay_vars().len()
            ),
            "int n_fixed_delays = n_delay_vars - n_inferred_delays;".into(),
            format!(
                "int n_transition_probabilities = {};", self.transitions_with_probabilities().len()
            ),
            format!(
                "int n_substates = {};",
                {
                    let mut pieces = Vec::new();
                    
                    for (state, isgamma) in self.state_isgamma_pairs() {
                        if isgamma {
                           pieces.push(format!("{}_gamma_shape", state));
                        }
                        else {
                            pieces.push("1".into());
                        }
                    }
                    
                    for obs_var in self.observation_variables() {
                        pieces.push(format!("{}_gamma_shape + 1", obs_var));
                    }
                    
                    pieces.join(" + ")
                }
            ),
            "int x_i[1 + n_delay_vars];".into(),
            "real x_r[1 + n_fixed_delays + n_transition_probabilities];".into(),
        ].join("\n")
    }
    
    fn assign_x_i(&self) -> String {
        let increment: String = "index += 1;".into();
        let mut sections = Vec::new();
    
        sections.push("int index = 1;".into());
    
        sections.push(vec![
            "x_i[index] = n_substates;".into(),
            increment.clone()
        ].join("\n"));
    
        for var in self.variables_with_delays(true, true) {
            sections.push(vec![
                format!("x_i[index] = {}_gamma_shape;", var),
                increment.clone(),
            ].join("\n"));
        }
    
        format_block(sections.join("\n\n"))
    }
    
    fn parameters_code(&self) -> String {
        vec![
            self.standard_declarations(false, true, true),
            self.mean_duration_declarations(false, true, true),
        ].join("\n\n")
    }
    
    fn transformed_parameters_code(&self) -> String {
        "".into()
    }
    
    fn model_code(&self) -> String {
        vec![
            "b ~ normal(0, 2);".into(),
            self.mean_duration_sampling_statements().join("\n"),
        ].join("\n\n")
    }
    
    fn mean_duration_sampling_statements(&self) -> Vec<String> {
        self.inferred_delay_vars().iter().map(|var| {
            format!("{}_mean_duration ~ normal(0, 20);", var)
        }).collect()
    }
    
    fn generated_quantities_code(&self) -> String {
        vec![
            self.gq_state_declarations(),
            self.gq_body(),
        ].join("\n\n")
    }
    
    fn gq_state_declarations(&self) -> String {
        let mut lines = Vec::new();
        
        for state in &self.structure.states {
            lines.push(format!("vector[n_times] {};", state.name));
        }
        
        for obs_var in self.observation_variables() {
            lines.push(format!("vector[n_times] c_{}_hidden;", obs_var));
        }
        
        lines.join("\n")
    }
    
    fn gq_body(&self) -> String {
        format_block(vec![
            vec![
                String::from("real initial_state[n_substates];"),
                "real params[1 + n_inferred_delays];".into(),
                "real ode_out[n_times, n_substates];".into(),
            ].join("\n"),
            self.extract_or_assign("params", self.params_entries(), false),
            self.gq_ics(),
            vec![
                String::from("ode_out = integrate_ode_rk45("),
                "  ode_ddt, initial_state, start_time, times, params, x_r, x_i".into(),
                ");".into()
            ].join("\n"),
            self.gq_assign_states(),
        ].join("\n\n"))
    }
    
    fn gq_ics(&self) -> String {
        let mut sections = Vec::new();
    
        sections.push("int index = 1;".into());
    
        for (state, isgamma) in self.state_isgamma_pairs() {
            if isgamma {
                sections.push(vec![
                    format!("initial_state[index:(index + {}_gamma_shape - 1)] = rep_array(", state),
                    format!("  {0}_init / {0}_gamma_shape, {0}_gamma_shape", state),
                    ");".into(),
                    format!("index += {}_gamma_shape;", state),
                ].join("\n"));
            }
            else {
                sections.push(vec![
                    format!("initial_state[index] = {}_init;", state),
                    "index += 1;".into(),
                ].join("\n"));
            }
        }
    
        for obs_var in self.observation_variables() {
            sections.push(vec![
                format!("initial_state[index:(index + {}_gamma_shape)] = rep_array(", obs_var),
                format!("  0.0, {}_gamma_shape + 1", obs_var),
                ");".into(),
                format!("index += {}_gamma_shape + 1;", obs_var),
            ].join("\n"));
        }
    
        format_block(sections.join("\n\n"))
    }
    
    fn gq_assign_states(&self) -> String {
        let mut sections = Vec::new();
        
        sections.push("int index = 1;".into());
        
        for (state, isgamma) in self.state_isgamma_pairs() {
            if isgamma {
                sections.push(vec![
                    format!(
                        "{0}[i] = sum(ode_out[i, index:(index + {0}_gamma_shape - 1)]);",
                        state
                    ),
                    format!("index += {}_gamma_shape;", state),
                ].join("\n"))
            }
            else {
                sections.push(vec![
                    format!("{0}[i] = ode_out[i, index];", state),
                    "index += 1;".into()
                ].join("\n"))
            }
        }
        
        for obs_var in self.observation_variables() {
            sections.push(vec![
                format!("c_{0}_hidden[i] = ode_out[i, index + {0}_gamma_shape];", obs_var),
                format!("index += {}_gamma_shape + 1;", obs_var),
            ].join("\n"))
        }
        
        sections.push({
            let mut lines = Vec::new();
        
            lines.push(format!("{}[i] = N", self.susceptible_state()));
        
            for (name, isgamma) in self.state_isgamma_pairs() {
                lines.push(format!("  - {}[i]", name));
            }
        
            lines.push(";".into());
        
            lines.join("\n")
        });
        
        
        format!("for(i in 1:n_times) {{\n{}\n}}", indent(&sections.join("\n\n"), 2, 1))
    }
}

fn indent(s: &str, n_spaces: usize, n_indent: usize) -> String {
    s.lines().map(|line| {
        format!("{}{}", (" ".repeat(n_spaces)).repeat(n_indent), line)
    }).collect::<Vec<_>>().join("\n")
}

fn generate_stan_code(
    functions: Vec<String>,
    data: String,
    transformed_data: String,
    parameters: String,
    transformed_parameters: String,
    model: String,
    generated_quantities: String,
) -> String {
    format!("{functions}\n\
        {data}\n\
        {transformed_data}\n\
        {parameters}\n\
        {transformed_parameters}\n\
        {model}\n\
        {generated_quantities}\n",
        functions = functions_code(functions),
        data = data_code(data),
        transformed_data = transformed_data_code(transformed_data),
        parameters = parameters_code(parameters),
        transformed_parameters = transformed_parameters_code(transformed_parameters),
        model = model_code(model),
        generated_quantities = generated_quantities_code(generated_quantities),
    )
}

fn format_block(body: String) -> String {
    format!("{{\n{}\n}}", indent(&body, 2, 1))
}

fn format_section(name: &str, body: String) -> String {
    format!(
        "{} {{\n{}\n}}\n",
        name,
        indent(&body, 2, 1)
    )
}

pub fn format_function(declaration: &str, body: &str) -> String {
    format!(
        "{declaration} {{\n\
            {body}\n\
            }}",
        declaration = declaration,
        body = indent(body, 2, 1)
    )
}

fn functions_code(functions: Vec<String>) -> String {
    format_section("functions", functions.join("\n\n"))
}

fn data_code(body: String) -> String {
    format_section("data", body)
}

fn transformed_data_code(body: String) -> String {
    format_section("transformed data", body)
}

fn parameters_code(body: String) -> String {
    format_section("parameters", body)
}

fn transformed_parameters_code(body: String) -> String {
    format_section("transformed parameters", body)
}

fn model_code(body: String) -> String {
    format_section("model", body)
}

fn generated_quantities_code(body: String) -> String {
    format_section("generated quantities", body)
}


#[derive(Debug, Serialize, Deserialize)]
pub enum Error {
    InvalidInputPath(String),
    InvalidInputFile(String),
    InputReadFailure,
    InvalidJson(JsonError),
}

#[derive(Debug, Serialize, Deserialize)]
pub struct JsonError {
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
pub struct InputData {
    structure: ModelStructure,
    config: ModelConfig,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct OutputData {
    stan_code: String,
}

pub fn run() -> Result<OutputData, Error> {
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
//    eprintln!("{}", stan_code);
    
    Ok(OutputData {
        stan_code
    })
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

#[cfg(test)]
mod test {
    fn run_test_file(name: &str) {
        use crate::stan::*;
        
        let input_filename = format!("tests/{}.json", name);
        let goal_filename = format!("tests/{}.stan", name);
        let gen_filename = format!("tests/{}-generated.stan", name);
    
        let json_data = read_data_from_file(&input_filename).unwrap();
        let input_data: InputData = serde_json::from_str(&json_data).unwrap();
        
        let goal_code = read_data_from_file(&goal_filename).unwrap();
        let gen_code = StanModel::new(
            input_data.structure,
            input_data.config
        ).generate_stan_code();
        write_data_to_file(&gen_filename, &gen_code);
        
        let mut goal_itr = goal_code.lines();
        let mut gen_itr = gen_code.lines();
        
        for (i, goal_line) in goal_itr.enumerate() {
            match gen_itr.next() {
                Some(gen_line) => {
                    if goal_line.ne(gen_line) {
                        eprintln!("line {} not equal:\n\nOriginal:\n{}\nGenerated:\n{}", i + 1, goal_line, gen_line);
                    }
                    assert_eq!(gen_line, goal_line);
                }
                None => {
                    eprintln!("line {} missing", i);
                    assert!(false);
                }
            };
        }
    }
    
    #[test]
    fn seir_with_delays() {
        run_test_file("seir-with-delays");
    }
}

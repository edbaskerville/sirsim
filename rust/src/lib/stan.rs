use crate::spec::{ModelStructure, ModelConfig};

pub struct StanModel {
}

impl StanModel {
    pub fn new(structure: ModelStructure, config: ModelConfig) -> Self {
        Self {}
    }
    
    pub fn generate_stan_code(&self) -> String {
        generate_stan_code(
            vec![],
            vec![],
            "".into(),
            vec![],
            "".into(),
            "".into(),
            "".into(),
        )
    }
}

fn indent(s: &str, n_spaces: usize, n_indent: usize) -> String {
    s.lines().map(|line| {
        format!("{}{}", (" ".repeat(n_spaces)).repeat(n_indent), line)
    }).collect::<Vec<_>>().join("\n")
}

fn generate_stan_code(
    functions: Vec<String>,
    data: Vec<String>,
    transformed_data: String,
    parameters: Vec<String>,
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

fn format_section(name: &str, body: String) -> String {
    format!(
        "{} {{\n\
            {}\n\
            }}\n\
        ",
        name,
        indent(&body, 2, 1)
    )
}

fn functions_code(functions: Vec<String>) -> String {
    format_section("functions", functions.join("\n\n"))
}

fn data_code(data_variables: Vec<String>) -> String {
    format_section("data", data_variables.join("\n"))
}

fn transformed_data_code(body: String) -> String {
    format_section("transformed data", body)
}

fn parameters_code(parameters: Vec<String>) -> String {
    format_section("parameters", parameters.join("\n"))
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

sirstan_generate <- function(
  sirstan_root,
  model_structure,
  config,
  sirstan_input_data_output_path = NULL,
  stan_data_output_path = NULL,
  stan_code_output_path = NULL
) {
  library(jsonlite)
  
  input_data_jsonobj <- sirstan_prepare_jsonobj_(model_structure, config)
  input_data_json <- sirstan_to_json(input_data_jsonobj)
  
  sirstan_exec_path <- file.path(sirstan_root, 'rust/target/debug/sirstan')
  
  output <- system(
    sirstan_exec_path,
    input = input_data_json,
    intern = TRUE
  )
  
  list(
    sirstan_output = output,
    sirstan_input_data_jsonobj = input_data_jsonobj,
    sirstan_input_data_json = input_data_json
  )
}

sirstan_to_json <- function(x) {
  library(jsonlite)
  toJSON(
    x,
    null = 'null',
    na = 'null',
    digits = NA,
    pretty = TRUE
  )
}

sirstan_prepare_jsonobj_ <- function(model_structure, config) {
  library(jsonlite)
  
  model_structure_ <- with(model_structure, list(
    susceptible_state = unbox(susceptible_state),
    states = {
      state_names <- names(states)
      lapply(1:length(states), function(i) {
        state <- states[[i]]
        list(
          name = unbox(state_names[i]),
          infectious = unbox(state$infectious),
          next_states = state$next_states
        )
      })
    },
    observation_variables = {
      ov_names <- names(observation_variables)
      lapply(1:length(observation_variables), function(i) {
        ov <- observation_variables[[i]]
        list(
          name = unbox(ov_names[i]),
          start_state = unbox(ov$start_state),
          end_state = unbox(ov$end_state)
        )
      })
    }
  ))
  
  format_parameter_distribution <- function(d) {
    if(is.character(d)) {
      list(
        tag = unbox('RawCode'),
        content = unbox(d)
      )
    }
    else {
      list(
        tag = unbox('Structured'),
        content = list(
          type_ = d$type,
          parameters = lapply(d$parameters, function(param) {
            format_parameter(param)
          })
        )
      )
    }
  }
  
  format_parameter <- function(param) {
    if(is.numeric(param)) {
      list(
        tag = unbox('Fixed'),
        content = unbox(param)
      )
    }
    else {
      list(
        tag = unbox('Inferred'),
        content = format_parameter_distribution(param)
      )
    }
  }
  
  format_delay <- function(delay_config) {
    with(delay_config, list(
      mean_duration = format_parameter(mean_duration),
      gamma_shape = format_parameter(gamma_shape)
    ))
  }
  
  format_observation_values <- function(values) {
    lapply(1:nrow(values), function(i) {
      c(values$start_time[i], values$end_time[i], values$value[i])
    })
  }
  
  format_observation_distribution <- function(d) {
    with(d, list(
      tag = unbox(type),
      content = lapply(parameters, format_parameter)
    ))
  }
  
  config_ <- with(config, {
    list(
      population_size = unbox(population_size),
      transmission_covariates = {
        tc_names <- names(transmission_covariates)
        lapply(1:length(transmission_covariates), function(i) {
          tc <- transmission_covariates[[i]]
          list(
            name = unbox(tc_names[i]),
            values = lapply(1:nrow(tc), function(j) {
              c(tc$time[j], tc$value[j])
            })
          )
        })
      },
      process_delays = lapply(process_delays, format_delay),
      observations = lapply(observations, function(obs) {
        with(obs, list(
          delay = format_delay(delay),
          distribution = format_observation_distribution(distribution),
          values = format_observation_values(values)
        ))
      })
    )
  })
  
  list(
    structure = model_structure_,
    config = config_
  )
}

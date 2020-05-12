sirsim_ibm_make_config <- function(
  rng_seed,
  output_path,
  write_to_stdout,
  record_all_events,
  
  t_final,
  
  n_ageclasses,
  susceptible_state,
  initial_infected_state,
  final_states,
  infected_states,
  contact_parameters,
  initial_counts,
  
  config_path = NULL
) {
  library(jsonlite)
  
  process_infected_state <- function(state) {
    with(state, list(
      name = unbox(name),
      infectious = unbox(infectious),
      mean_duration = unbox(mean_duration),
      gamma_shape = unbox(gamma_shape),
      next_states = next_states,
      probabilities = state$probabilities
    ))
  }
  
  process_contact_parameters_item <- function(cp_item) {
    list(
      beta = unbox(cp_item$beta),
      C = cp_item$C,
      t_end = unbox(cp_item$t_end)
    )
  }
  
  config <- list(
    rng_seed = unbox(rng_seed),
    output_path = unbox(output_path),
    write_to_stdout = unbox(write_to_stdout),
    record_all_events = unbox(record_all_events),
    t_final = unbox(t_final),
    
    n_ageclasses = unbox(n_ageclasses),
    susceptible_state = unbox(susceptible_state),
    initial_infected_state = unbox(initial_infected_state),
    final_states = final_states,
    
    infected_states = lapply(infected_states, process_infected_state),
    contact_parameters = lapply(contact_parameters, process_contact_parameters_item),
    initial_counts = initial_counts
  )
  
  config_json <- toJSON(
    config, null = 'null', pretty = TRUE,
    matrix = 'rowmajor',
    digits = NA
  )
  
  if(!is.null(config_path)) {
    write(config_json, config_path)
  }
  
  config
}

sirsim_ibm_simulate <- function(
  sirsim_root,
  config,
  with_rust_optimizations = TRUE,
  with_rust_backtrace = FALSE
) {
  library(jsonlite)
  
  sirsim_build(sirsim_root)
  
  exec_path <- normalizePath(
    file.path(
      sirsim_root, 'rust', 'target',
      if(with_rust_optimizations) 'release' else 'debug',
      'sirsim'
    )
  )
  
  backtrace <- if(with_rust_backtrace) 'RUST_BACKTRACE=full ' else ''
  config_path <- if(is.list(config)) {
    ''
  } else {
    sprintf('"%s"', normalizePath(config))
  }
  
  input <- if(is.list(config)) {
    toJSON(
      config, null = 'null', pretty = TRUE,
      matrix = 'rowmajor',
      digits = NA
    )
  } else NULL
  
  output <- system(sprintf(
    '%s"%s" %s',
    backtrace,
    exec_path,
    config_path
  ), intern = TRUE, input = input)
  
  output
}

sirsim_build <- function(
  sirsim_root,
  with_rust_optimizations = TRUE
) {
  old_wd <- normalizePath(getwd())
  rust_root <- normalizePath(file.path(sirsim_root, 'rust'))
  
  setwd(rust_root)
  err <- system(sprintf(
    '~/.cargo/bin/cargo build%s',
    if(with_rust_optimizations) ' --release' else ''
  ))
  setwd(old_wd)
  
  if(err != 0) {
    stop('Failed to build sirsim.')
  }
}

sirsim_compute_Rt <- function(con) {
  library(dplyr)
  
  df_with_gaps <- sirsim_get_infection_counts_by_id(con) %>%
    collect() %>%
    mutate(t_discrete = ceiling(t_inf_source)) %>%
    group_by(t_discrete) %>%
    summarize(n_infected = sum(n_infected), n = n(), Rt = n_infected / n) %>%
    ungroup() %>%
    select(t = t_discrete, n_infected, n, Rt)
  
  t_start <- ceiling(min(df_with_gaps$t))
  t_end <- ceiling(max(df_with_gaps$t))
  
  tibble(t = t_start:t_end) %>%
    left_join(
      df_with_gaps, by = 't'
    )
}

sirsim_get_infection_counts_by_id <- function(con) {
  library(dplyr)
  
  inf_tbl <- tbl(con, 'Infections')
  
  all_sources <- inf_tbl %>%
    select(t_inf_source = t, source_id = infected_id)
  all_infecteds <- inf_tbl %>%
    select(source_id = infectious_id, infected_id = infected_id)
  
  infection_counts <- all_sources %>%
    inner_join(
      all_infecteds,
      by = 'source_id'
    ) %>%
    group_by(source_id, t_inf_source) %>%
    summarize(n_infected = n()) %>%
    ungroup()
  
  result <- all_sources %>%
    left_join(
      infection_counts %>% select(source_id, n_infected),
      by = c('source_id')
    ) %>%
    mutate(n_infected = ifelse(is.na(n_infected), 0, n_infected))
  result
}

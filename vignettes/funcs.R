# Common functions for analysis
train_mofa <- function(mae, num_factors, convergence_mode) {
  # Create MOFA model
  model <- create_mofa_from_MultiAssayExperiment(
    mae,
    extract_metadata = TRUE
  )

  # Set model's options
  model_opts <- get_default_model_options(model)
  model_opts$num_factors <- num_factors
  train_opts <- get_default_data_options(model)
  # train_opts$convergence_mode <- convergence_mode

  # Prepare MOFA model
  model <- prepare_mofa(
    object = model,
    model_options = model_opts,
    # training_options = train_opts
  )

  # Train model
  model <- run_mofa(model, use_basilisk = TRUE)
  
  return(model)
}
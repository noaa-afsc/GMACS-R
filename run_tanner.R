script_path <- function(default) {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
  }
  source_file <- sys.frames()[[1]]$ofile
  if (!is.null(source_file)) {
    return(normalizePath(source_file, mustWork = TRUE))
  }
  normalizePath(default, mustWork = TRUE)
}

repo_root <- dirname(script_path("run_tanner.R"))
source(file.path(repo_root, "R", "read_gmacs_bbrkc.R"))
source(file.path(repo_root, "R", "read_gmacs_tanner.R"))

input_root <- Sys.getenv("GMACS_TANNER_ROOT", file.path(repo_root, "examples", "tanners"))
inputs <- read_tanner_inputs(input_root)
validate_tanner_inputs(inputs)

cat("Tanner GMACS inputs parsed\n")
cat("  root:", input_root, "\n")
cat("  years:", inputs$data$dimensions[["syr"]], "-", inputs$data$dimensions[["nyr"]], "\n")
cat("  fleets:", inputs$data$dimensions[["nfleet"]], "\n")
cat("  size classes:", inputs$data$dimensions[["nclass"]], "\n")
cat("  fishing fleets:", paste(inputs$data$fishing_fleet_names, collapse = ", "), "\n")
cat("  surveys:", paste(inputs$data$survey_names, collapse = ", "), "\n")
cat("  catch dataframes:", inputs$data$catch$n, " rows:", sum(inputs$data$catch$rows), "\n")
cat("  survey dataframes:", inputs$data$survey$n, " rows:", sum(inputs$data$survey$rows), "\n")
cat("  size-composition dataframes:", inputs$data$size_comp$n, " rows:", sum(inputs$data$size_comp$rows), "\n")
cat("  growth observations:", inputs$data$growth$n, "\n")
cat("  maturity observations:", inputs$data$maturity$n, "\n")
cat("  parameter rows:", nrow(inputs$parameters), "\n")
cat("  active estimated parameters:", sum(inputs$parameters$active), "\n")
if (!inputs$admb_reference_available) {
  cat("  ADMB parity outputs: not available in this folder; likelihood audit not run.\n")
}

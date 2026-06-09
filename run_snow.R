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

repo_root <- dirname(script_path("run_snow.R"))
source(file.path(repo_root, "R", "read_gmacs_bbrkc.R"))
source(file.path(repo_root, "R", "gmacs_rtmb_nll.R"))
source(file.path(repo_root, "R", "read_gmacs_snow.R"))

input_root <- Sys.getenv("GMACS_SNOW_ROOT", file.path(repo_root, "examples", "snow"))
inputs <- read_snow_inputs(input_root)
validate_snow_inputs(inputs)

cat("Snow crab GMACS inputs parsed\n")
cat("  root:", input_root, "\n")
cat("  years:", inputs$data$dimensions[["syr"]], "-", inputs$data$dimensions[["nyr"]], "\n")
cat("  fleets:", inputs$data$dimensions[["nfleet"]], "\n")
cat("  size classes:", inputs$data$dimensions[["nclass"]], "\n")
cat("  catch rows:", sum(inputs$data$catch$rows), "\n")
cat("  survey rows:", nrow(inputs$data$survey$data), "\n")
cat("  size-composition rows:", sum(inputs$data$size_comp$rows), "\n")
cat("  parameter blocks:", paste(names(inputs$parameters), collapse = ", "), "\n")
cat("  ADMB likelihood entries:", length(inputs$admb_reference$nloglike), "\n")
cat("  ADMB penalty entries:", length(inputs$admb_reference$nlog_penalty), "\n")
cat("  ADMB prior entries:", length(inputs$admb_reference$prior_density), "\n")
cat("  ADMB catch-fit rows:", nrow(inputs$admb_catch_fit), "\n")
cat("  ADMB index-fit rows:", nrow(inputs$admb_index_fit), "\n")
cat("  ADMB size-fit rows:", nrow(inputs$admb_size_fit), "\n")

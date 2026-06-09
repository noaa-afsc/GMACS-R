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

repo_root <- dirname(script_path("run_bbrkc.R"))
source(file.path(repo_root, "R", "read_gmacs_bbrkc.R"))
source(file.path(repo_root, "R", "gmacs_rtmb_nll.R"))

input_root <- Sys.getenv("GMACS_BBRKC_ROOT", file.path(repo_root, "build", "BBRKC"))
required_inputs <- file.path(input_root, c("gmacs.dat", "gmacs.pin", "Gmacsall.out", "gmacs.rep"))
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop(
    "Missing BBRKC input/reference files under ", input_root, ":\n  ",
    paste(basename(missing_inputs), collapse = "\n  "),
    "\nSet GMACS_BBRKC_ROOT to the directory containing these files.",
    call. = FALSE
  )
}

inputs <- read_bbrkc_inputs(input_root)
validate_gmacs_inputs(inputs)

cat("BBRKC GMACS inputs parsed\n")
cat("  years:", inputs$data$dimensions[["syr"]], "-", inputs$data$dimensions[["nyr"]], "\n")
cat("  fleets:", inputs$data$dimensions[["nfleet"]], "\n")
cat("  size classes:", inputs$data$dimensions[["nclass"]], "\n")
cat("  catch rows:", sum(inputs$data$catch$rows), "\n")
cat("  survey rows:", nrow(inputs$data$survey$data), "\n")
cat("  size-composition rows:", sum(inputs$data$size_comp$rows), "\n")
cat("  parameter blocks:", paste(names(inputs$parameters), collapse = ", "), "\n")

if (requireNamespace("RTMB", quietly = TRUE)) {
  map <- make_bbrkc_rtmb_map(inputs)
  cat("  mapped active parameters:", count_active_map_parameters(map), "\n")
  obj <- make_gmacs_rtmb_object(inputs)
  cat("RTMB object taped\n")
  cat("  active parameter length:", length(obj$par), "\n")
  cat("  initial objective:", obj$fn(obj$par), "\n")
  report <- obj$report()
  audit <- audit_bbrkc_admb_components(report, inputs)
  cat("  ADMB likelihood/penalty/prior audit:\n")
  print_bbrkc_admb_audit(audit)
  cat("  weighted penalty sum:", signif(sum(report$weighted_nlog_penalty), 8), "\n", sep = "")
  model <- initialize_bbrkc_model_parameters(make_gmacs_parameter_list(inputs), make_gmacs_rtmb_data(inputs))
  for (component in c("capture", "retained", "discard")) {
    cmp <- compare_selectivity_to_admb(model, inputs, component)
    cat("  ", component, " selectivity max abs diff:", signif(cmp$max_abs_diff, 4), "\n", sep = "")
  }
  f_cmp <- compare_fully_selected_f_to_admb(model, inputs)
  cat("  fully selected F max abs diff:", signif(f_cmp$max_abs_diff, 4), "\n", sep = "")
  f_size_cmp <- compare_f_at_size_to_admb(model, inputs)
  cat("  F at size max abs diff:", signif(f_size_cmp$max_abs_diff, 4), "\n", sep = "")
  m_cmp <- compare_natural_mortality_to_admb(model, inputs)
  cat("  natural mortality max abs diff:", signif(m_cmp$max_abs_diff, 4), "\n", sep = "")
  z_cmp <- compare_total_mortality_to_admb(model, inputs, "continuous")
  cat("  continuous total mortality max abs diff:", signif(z_cmp$max_abs_diff, 4), "\n", sep = "")
  z2_cmp <- compare_total_mortality_to_admb(model, inputs, "discrete")
  cat("  discrete total mortality max abs diff:", signif(z2_cmp$max_abs_diff, 4), "\n", sep = "")
  molt_cmp <- compare_molt_probability_to_admb(model, inputs)
  cat("  molt probability max abs diff:", signif(molt_cmp$max_abs_diff, 4), "\n", sep = "")
  growth_cmp <- compare_growth_transition_to_admb(model, inputs)
  cat("  growth transition max abs diff:", signif(growth_cmp$max_abs_diff, 4), "\n", sep = "")
  recruit_cmp <- compare_recruits_to_admb(model, inputs)
  cat("  recruits max abs diff:", signif(recruit_cmp$max_abs_diff, 4), "\n", sep = "")
  for (component in c(
    "total", "males", "females", "males_new", "females_new",
    "males_old", "females_old"
  )) {
    n_cmp <- compare_numbers_at_size_to_admb(model, inputs, component)
    cat(
      "  N(", component, ") max abs diff:",
      signif(n_cmp$max_abs_diff, 4), "\n",
      sep = ""
    )
  }
  catch_cmp <- compare_catch_fit_to_admb(model, inputs)
  cat("  catch predicted max abs diff:", signif(catch_cmp$predicted_max_abs_diff, 4), "\n", sep = "")
  cat("  catch residual max abs diff:", signif(catch_cmp$residual_max_abs_diff, 4), "\n", sep = "")
  index_cmp <- compare_index_fit_to_admb(model, inputs)
  cat("  index predicted max abs diff:", signif(index_cmp$predicted_max_abs_diff, 4), "\n", sep = "")
  cat("  index q max abs diff:", signif(index_cmp$q_max_abs_diff, 4), "\n", sep = "")
  size_cmp <- compare_size_composition_to_admb(model, inputs)
  cat("  size composition predicted max abs diff:", signif(size_cmp$max_abs_diff, 4), "\n", sep = "")
} else {
  cat("RTMB is not installed; parser and parameter construction were checked only.\n")
}

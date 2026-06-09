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

repo_root <- dirname(script_path("optimize_bbrkc.R"))
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
obj <- make_gmacs_rtmb_object(inputs)

initial_objective <- obj$fn(obj$par)
initial_gradient <- obj$gr(obj$par)

cat("BBRKC RTMB optimization smoke test\n")
cat("  active parameter length:", length(obj$par), "\n")
cat("  initial objective:", initial_objective, "\n")
cat("  finite gradient:", all(is.finite(initial_gradient)), "\n")
cat("  max abs gradient:", signif(max(abs(initial_gradient)), 6), "\n")

opt <- nlminb(
  start = obj$par,
  objective = obj$fn,
  gradient = obj$gr,
  control = list(eval.max = 10, iter.max = 5)
)

cat("  optimizer convergence:", opt$convergence, "\n")
cat("  optimizer message:", opt$message, "\n")
cat("  final objective:", opt$objective, "\n")
cat("  objective change:", opt$objective - initial_objective, "\n")

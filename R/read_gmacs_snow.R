read_snow_inputs <- function(root = "examples/snow") {
  inputs <- read_bbrkc_inputs(root)
  inputs$example_name <- "snow"
  inputs
}

validate_snow_inputs <- function(inputs) {
  validate_gmacs_inputs(inputs)
}

snow_likelihood_summary <- function(inputs) {
  data.frame(
    block = c("nloglike", "nlogPenalty", "priorDensity"),
    entries = c(
      length(inputs$admb_reference$nloglike),
      length(inputs$admb_reference$nlog_penalty),
      length(inputs$admb_reference$prior_density)
    ),
    sum = c(
      sum(inputs$admb_reference$nloglike),
      sum(inputs$admb_reference$nlog_penalty),
      sum(inputs$admb_reference$prior_density)
    ),
    stringsAsFactors = FALSE
  )
}

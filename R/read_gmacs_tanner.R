tanner_file_set <- function(root = "examples/tanners") {
  files <- list(
    main = "gmacs_26_22_03d5_aEffExp1.Run1.dat",
    data = "TannerCrab_Data202509.EffExp1.dat",
    control = "TannerCrab_26_22_03d5_a.Run1.ctl",
    project = "TannerCrab_26_22_03d5_a.prj",
    par = "gmacs_26_22_03d5_aEffExp1.par"
  )
  lapply(files, function(file) file.path(root, file))
}

tanner_admb_report_candidates <- function(root = "examples/tanners") {
  list(
    all = file.path(root, c(
      "Gmacsall.out",
      "gmacsall.out",
      "GmacsAll.out",
      "TannerCrab_Gmacsall.out",
      "gmacs_26_22_03d5_aEffExp1.Gmacsall.out",
      "gmacs_26_22_03d5_aEffExp1.all.out"
    )),
    report = file.path(root, c(
      "gmacs.rep",
      "Gmacs.rep",
      "TannerCrab.rep",
      "TannerCrab_26_22_03d5_a.rep",
      "gmacs_26_22_03d5_aEffExp1.rep"
    ))
  )
}

tanner_admb_report_files <- function(root = "examples/tanners") {
  candidates <- tanner_admb_report_candidates(root)
  files <- lapply(candidates, function(x) {
    hit <- x[file.exists(x)]
    if (length(hit)) hit[1] else NA_character_
  })
  names(files) <- names(candidates)
  files
}

read_tanner_admb_reference <- function(root = "examples/tanners") {
  files <- tanner_admb_report_files(root)
  out <- list(files = files)

  if (!is.na(files$all)) {
    out$likelihood <- read_admb_likelihood_reference(files$all)
    out$catch_fit <- read_admb_catch_fit_summary(files$all)
    out$index_fit <- read_admb_index_fit_summary(files$all)
    out$size_fit <- read_admb_size_fit_summary(files$all)
    out$overall_summary <- read_admb_overall_summary(files$all)
  }

  if (!is.na(files$report)) {
    out$selectivity <- list(
      capture = read_admb_selectivity_block(files$report, "slx_capture"),
      retained = read_admb_selectivity_block(files$report, "slx_retaind"),
      discard = read_admb_selectivity_block(files$report, "slx_discard")
    )
  }

  out$available <- !is.na(files$all) && !is.na(files$report)
  out
}

first_numeric_on_line <- function(lines, pattern) {
  hit <- grep(pattern, lines)
  if (!length(hit)) {
    stop("Could not find Tanner label: ", pattern, call. = FALSE)
  }
  numeric_values(lines[hit[1]])[1]
}

next_numeric_line <- function(lines, index) {
  j <- index + 1L
  while (j <= length(lines)) {
    values <- numeric_values(lines[j])
    if (length(values)) {
      return(values)
    }
    j <- j + 1L
  }
  numeric()
}

read_tanner_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  values <- strip_comment(lines)
  values <- values[nzchar(values)]
  values <- values[!grepl("^[=-]+\\s*$", values)]
  jitter <- numeric_values(values[7])

  list(
    data_file = values[1],
    control_file = values[2],
    project_file = values[3],
    weight_unit = values[4],
    numbers_unit = values[5],
    stock_name = values[6],
    is_jittered = as.integer(jitter[1]),
    use_pin_flag = as.integer(jitter[2]),
    sd_jitter = jitter[3],
    output_variance = as.integer(vapply(values[8:12], numeric_values, numeric(1))),
    nyr_retro = as.integer(numeric_values(values[13])),
    turn_off_phase = as.integer(numeric_values(values[14])),
    stop_after_fn_call = as.integer(numeric_values(values[15])),
    calc_ref_points = as.integer(numeric_values(values[16])),
    use_pin_file = as.integer(numeric_values(values[17])),
    verbose = as.integer(numeric_values(values[18]))
  )
}

rows_between <- function(lines, start_pattern, end_pattern = NULL) {
  start <- grep(start_pattern, lines)
  if (!length(start)) {
    return(integer())
  }
  if (is.null(end_pattern)) {
    section <- lines[start[1]:length(lines)]
  } else {
    end <- grep(end_pattern, lines)
    end <- end[end > start[1]]
    if (!length(end)) {
      section <- lines[start[1]:length(lines)]
    } else {
      section <- lines[start[1]:(end[1] - 1L)]
    }
  }
  as.integer(vapply(grep("#--number of rows in dataframe", section, value = TRUE), function(x) {
    numeric_values(x)[1]
  }, numeric(1)))
}

section_between <- function(lines, start_pattern, end_pattern = NULL) {
  start <- grep(start_pattern, lines)
  if (!length(start)) {
    return(character())
  }
  if (is.null(end_pattern)) {
    return(lines[start[1]:length(lines)])
  }
  end <- grep(end_pattern, lines)
  end <- end[end > start[1]]
  if (!length(end)) {
    return(lines[start[1]:length(lines)])
  }
  lines[start[1]:(end[1] - 1L)]
}

read_tanner_survey_observations <- function(lines) {
  section <- section_between(lines, "RELATIVE ABUNDANCE DATA", "SIZE COMPS")
  starts <- grep("^##--[0-9]+ fleet:", section)
  records <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    pos <- starts[i]
    metadata <- strip_comment(section[(pos + 1L):(pos + 6L)])
    n_rows <- as.integer(numeric_values(section[pos + 7L])[1])
    data_start <- pos + 9L
    rows <- numeric_lines(section[data_start:(data_start + n_rows - 1L)])
    frame <- as.data.frame(as_matrix_rows(rows, 7, paste0("tanner_survey_", i)))
    names(frame) <- c("q_index", "year", "season", "obs", "cv", "multiplier", "cpue_time")
    frame$series <- i
    frame$units_type <- metadata[1]
    frame$index_type <- metadata[2]
    frame$fleet <- metadata[3]
    frame$sex <- metadata[4]
    frame$maturity <- metadata[5]
    frame$shell <- metadata[6]
    records[[i]] <- frame
  }
  do.call(rbind, records)
}

read_tanner_size_composition_observations <- function(lines) {
  section <- section_between(lines, "SIZE COMPS", "GROWTH DATA")
  starts <- grep("^##--fleet:", section)
  records <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    pos <- starts[i]
    metadata <- strip_comment(section[(pos + 1L):(pos + 5L)])
    n_rows <- as.integer(numeric_values(section[pos + 6L])[1])
    n_bins <- as.integer(numeric_values(section[pos + 7L])[1])
    data_start <- pos + 9L
    rows <- numeric_lines(section[data_start:(data_start + n_rows - 1L)])
    frame <- as.data.frame(as_matrix_rows(rows, 3 + n_bins, paste0("tanner_size_comp_", i)))
    names(frame) <- c("year", "season", "stage1_effn", paste0("bin_", seq_len(n_bins)))
    frame$series <- i
    frame$catch_type <- metadata[1]
    frame$fleet <- metadata[2]
    frame$sex <- metadata[3]
    frame$maturity <- metadata[4]
    frame$shell <- metadata[5]
    records[[i]] <- frame
  }
  do.call(rbind, records)
}

read_tanner_data_summary <- function(path) {
  lines <- readLines(path, warn = FALSE)
  size_breaks <- next_numeric_line(lines, grep("size_breaks", lines)[1])

  fleet_start <- grep("# Fishing fleet", lines)[1] + 1L
  survey_start <- grep("# Survey names", lines)[1]
  fleet_names <- strip_comment(lines[fleet_start:(survey_start - 1L)])
  fleet_names <- fleet_names[nzchar(fleet_names)]

  season_start <- grep("# Are the seasons", lines)[1]
  survey_names <- strip_comment(lines[(survey_start + 1L):(season_start - 1L)])
  survey_names <- survey_names[nzchar(survey_names)]

  growth_type <- first_numeric_on_line(lines, "#--growth type")
  growth_n <- first_numeric_on_line(lines, "#--number of observations \\(nobs_growth\\)")
  maturity_start <- grep("MALE MATURITY OGIVE DATA", lines)[1]
  maturity_n <- first_numeric_on_line(lines[maturity_start:length(lines)], "number of observations")

  list(
    dimensions = c(
      syr = as.integer(first_numeric_on_line(lines, "# Start")),
      nyr = as.integer(first_numeric_on_line(lines, "# End")),
      nseason = as.integer(first_numeric_on_line(lines, "Number of seasons")),
      nfleet = as.integer(first_numeric_on_line(lines, "Number of fleets")),
      nsex = as.integer(first_numeric_on_line(lines, "Number of sexes")),
      nshell = as.integer(first_numeric_on_line(lines, "Number of shell")),
      nmature = as.integer(first_numeric_on_line(lines, "Number of maturity")),
      nclass = as.integer(first_numeric_on_line(lines, "Number of size-classes"))
    ),
    season_recruitment = as.integer(first_numeric_on_line(lines, "Season recruitment")),
    season_growth = as.integer(first_numeric_on_line(lines, "Season molting")),
    season_ssb = as.integer(first_numeric_on_line(lines, "Season to calculate SSB")),
    season_N = as.integer(first_numeric_on_line(lines, "Season for N output")),
    n_size_sex = as.integer(c(
      next_numeric_line(lines, grep("# maximum size-class", lines)[1])[1],
      next_numeric_line(lines, grep("# maximum size-class", lines)[1] + 1L)[1]
    )),
    size_breaks = size_breaks,
    mid_points = size_breaks[-length(size_breaks)] + diff(size_breaks) / 2,
    m_prop_type = as.integer(next_numeric_line(lines, grep("Natural.*mortality.*input.*type", lines)[1])[1]),
    m_prop = as.numeric(vapply(grep("#--Season", lines, value = TRUE)[1:6], function(x) {
      numeric_values(x)[1]
    }, numeric(1))),
    fleet_names = c(fleet_names, survey_names),
    fishing_fleet_names = fleet_names,
    survey_names = survey_names,
    season_type = as.integer(next_numeric_line(lines, season_start)),
    catch = list(
      n = as.integer(first_numeric_on_line(lines, "#--number of catch dataframes")),
      rows = rows_between(lines, "CATCH DATA", "RELATIVE ABUNDANCE DATA")
    ),
    survey = list(
      n = as.integer(first_numeric_on_line(lines, "#--number of dataframes")),
      rows = rows_between(lines, "RELATIVE ABUNDANCE DATA", "SIZE COMPS")
    ),
    size_comp = list(
      n = as.integer(first_numeric_on_line(lines, "#--number of size comps dataframes")),
      rows = rows_between(lines, "SIZE COMPS", "GROWTH DATA")
    ),
    survey_observations = read_tanner_survey_observations(lines),
    size_comp_observations = read_tanner_size_composition_observations(lines),
    growth = list(type = as.integer(growth_type), n = as.integer(growth_n)),
    maturity = list(n = as.integer(maturity_n))
  )
}

read_tanner_parameter_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  records <- lapply(lines, function(line) {
    values <- numeric_values(line)
    if (!length(values) || !grepl("#", line, fixed = TRUE)) {
      return(NULL)
    }
    comment <- trimws(sub("^.*#", "", line))
    fields <- strsplit(comment, "\\s+")[[1]]
    if (!length(fields) || is.na(suppressWarnings(as.integer(fields[1])))) {
      return(NULL)
    }
    active <- !grepl("not_estd", comment, fixed = TRUE)
    est_index <- if (active) {
      tail(suppressWarnings(as.integer(fields)), 1)
    } else {
      NA_integer_
    }
    data.frame(
      value = values[1],
      parameter_index = as.integer(fields[1]),
      label = paste(fields[-1], collapse = " "),
      active = active,
      estimated_index = est_index,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, records[!vapply(records, is.null, logical(1))])
}

read_tanner_inputs <- function(root = "examples/tanners") {
  files <- tanner_file_set(root)
  missing <- unlist(files)[!file.exists(unlist(files))]
  if (length(missing)) {
    stop("Missing Tanner input files:\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
  }

  list(
    files = files,
    main = read_tanner_main(files$main),
    data = read_tanner_data_summary(files$data),
    parameters = read_tanner_parameter_table(files$par),
    admb_reference = read_tanner_admb_reference(root)
  )
}

validate_tanner_inputs <- function(inputs) {
  checks <- c(
    size_breaks = length(inputs$data$size_breaks) == inputs$data$dimensions[["nclass"]] + 1L,
    n_size_sex = length(inputs$data$n_size_sex) == inputs$data$dimensions[["nsex"]],
    catch_rows = length(inputs$data$catch$rows) == inputs$data$catch$n,
    survey_rows = length(inputs$data$survey$rows) == inputs$data$survey$n,
    size_comp_rows = length(inputs$data$size_comp$rows) == inputs$data$size_comp$n,
    growth_rows = inputs$data$growth$n > 0L,
    parameter_rows = nrow(inputs$parameters) > 0L
  )
  if (!all(checks)) {
    stop("Tanner validation failed: ", paste(names(checks)[!checks], collapse = ", "), call. = FALSE)
  }
  invisible(checks)
}

tanner_admb_reference_available <- function(inputs) {
  isTRUE(inputs$admb_reference$available)
}

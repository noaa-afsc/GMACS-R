strip_comment <- function(x) {
  trimws(sub("#.*$", "", x))
}

numeric_values <- function(x) {
  x <- strip_comment(x)
  if (!nzchar(x)) {
    return(numeric())
  }
  out <- suppressWarnings(scan(text = x, what = numeric(), quiet = TRUE))
  out
}

numeric_lines <- function(lines) {
  vals <- lapply(lines, numeric_values)
  vals[lengths(vals) > 0]
}

line_after <- function(lines, pattern) {
  i <- grep(pattern, lines, fixed = TRUE)
  if (!length(i)) {
    stop("Could not find label: ", pattern, call. = FALSE)
  }
  j <- i[1] + 1L
  while (j <= length(lines) && !nzchar(strip_comment(lines[j]))) {
    j <- j + 1L
  }
  lines[j]
}

section_lines <- function(lines, start_pattern, end_pattern) {
  start <- grep(start_pattern, lines, fixed = TRUE)
  end <- grep(end_pattern, lines, fixed = TRUE)
  if (!length(start)) {
    stop("Could not find section: ", start_pattern, call. = FALSE)
  }
  if (!length(end)) {
    stop("Could not find section end: ", end_pattern, call. = FALSE)
  }
  lines[(start[1] + 1L):(end[1] - 1L)]
}

as_matrix_rows <- function(rows, ncol, what) {
  bad <- lengths(rows) != ncol
  if (any(bad)) {
    stop(
      what, " expected ", ncol, " columns, but got ",
      paste(unique(lengths(rows)[bad]), collapse = ", "),
      call. = FALSE
    )
  }
  do.call(rbind, rows)
}

read_gmacs_main <- function(path) {
  lines <- readLines(path, warn = FALSE)
  values <- strip_comment(lines)
  values <- values[nzchar(values)]
  values <- values[!grepl("^[=-]+\\s*$", values)]

  list(
    data_file = values[1],
    control_file = values[2],
    project_file = values[3],
    weight_unit = values[4],
    numbers_unit = values[5],
    stock_name = values[6],
    is_jittered = as.integer(values[7]),
    sd_jitter = as.numeric(values[8]),
    output_variance = as.integer(vapply(values[9:13], numeric_values, numeric(1))),
    nyr_retro = as.integer(numeric_values(values[14])),
    turn_off_phase = as.integer(numeric_values(values[15])),
    stop_after_fn_call = as.integer(numeric_values(values[16])),
    calc_ref_points = as.integer(numeric_values(values[17])),
    use_pin_file = as.integer(numeric_values(values[18])),
    verbose = as.integer(numeric_values(values[19]))
  )
}

read_gmacs_data <- function(path) {
  lines <- readLines(path, warn = FALSE)

  dims <- c(
    syr = numeric_values(line_after(lines, "#_Start year")),
    nyr = numeric_values(line_after(lines, "#_End year")),
    nseason = numeric_values(line_after(lines, "#_Number of seasons")),
    nfleet = numeric_values(line_after(lines, "#_Number of fleets")),
    nsex = numeric_values(line_after(lines, "#_Number of sexes")),
    nshell = numeric_values(line_after(lines, "#_Number of shell condition types")),
    nmature = numeric_values(line_after(lines, "#_Number of maturity types")),
    nclass = numeric_values(line_after(lines, "#_Number of size-classes"))
  )
  nclass <- as.integer(dims[["nclass"]])
  nseason <- as.integer(dims[["nseason"]])
  nsex <- as.integer(dims[["nsex"]])
  nyr <- as.integer(dims[["nyr"]])
  syr <- as.integer(dims[["syr"]])

  n_size_sex <- numeric_values(line_after(lines, "#_maximum size-class"))
  size_breaks <- numeric_values(line_after(lines, "#_size_breaks"))
  season_recruitment <- as.integer(numeric_values(line_after(lines, "#_Season recruitment occurs"))[1])
  season_growth <- as.integer(numeric_values(line_after(lines, "#_Season molting and growth occurs"))[1])
  season_ssb <- as.integer(numeric_values(line_after(lines, "#_Season to calculate SSB"))[1])
  season_N <- as.integer(numeric_values(line_after(lines, "#_Season for N output"))[1])

  m_section <- section_lines(lines, "##_Natural mortality", "##_Fishery and survey definition")
  m_nums <- numeric_lines(m_section)
  m_prop_type <- as.integer(m_nums[[1]][1])
  m_prop_rows <- if (m_prop_type == 1L) 1L else nyr - syr + 1L
  m_prop <- as_matrix_rows(m_nums[2:(1 + m_prop_rows)], nseason, "m_prop")

  fishery_section <- section_lines(lines, "##_Fishery and survey definition", "##_Catch data")
  fishery_clean <- strip_comment(fishery_section)
  fishery_clean <- fishery_clean[nzchar(fishery_clean)]
  fleet_names <- strsplit(fishery_clean[1], "\\s+")[[1]]
  survey_names <- strsplit(fishery_clean[2], "\\s+")[[1]]
  season_type <- numeric_values(fishery_clean[3])

  catch_section <- section_lines(lines, "##_Catch data", "##_Relative abundance data")
  catch_nums <- numeric_lines(catch_section)
  catch_fmt <- as.integer(catch_nums[[1]][1])
  n_catch_df <- as.integer(catch_nums[[2]][1])
  n_catch_rows <- as.integer(catch_nums[[3]])
  catch_data <- as_matrix_rows(
    catch_nums[4:(3 + sum(n_catch_rows))],
    11,
    "catch_data"
  )
  catch_frames <- split(
    as.data.frame(catch_data),
    rep(seq_len(n_catch_df), n_catch_rows)
  )
  catch_frames <- lapply(catch_frames, `names<-`, c(
    "year", "season", "fleet", "sex", "obs", "cv", "type", "units",
    "mult", "effort", "discard_mortality"
  ))

  survey_section <- section_lines(lines, "##_Relative abundance data", "##_Size composition")
  survey_nums <- numeric_lines(survey_section)
  survey_fmt <- as.integer(survey_nums[[1]][1])
  n_surveys <- as.integer(survey_nums[[2]][1])
  index_type <- as.integer(survey_nums[[3]])
  n_survey_rows <- as.integer(survey_nums[[4]][1])
  survey_data <- as.data.frame(as_matrix_rows(
    survey_nums[5:(4 + n_survey_rows)],
    10,
    "survey_data"
  ))
  names(survey_data) <- c(
    "index", "year", "season", "fleet", "sex", "maturity", "obs",
    "cv", "units", "cpue_time"
  )

  comp_section <- section_lines(lines, "##_Size composition", "##_Growth data")
  comp_nums <- numeric_lines(comp_section)
  comp_fmt <- as.integer(comp_nums[[1]][1])
  n_size_comps <- as.integer(comp_nums[[2]][1])
  n_size_comp_rows <- as.integer(comp_nums[[3]])
  n_size_comp_cols <- as.integer(comp_nums[[4]])
  comp_rows <- comp_nums[5:(4 + sum(n_size_comp_rows))]
  comp_frames <- vector("list", n_size_comps)
  pos <- 1L
  for (i in seq_len(n_size_comps)) {
    rows <- comp_rows[pos:(pos + n_size_comp_rows[i] - 1L)]
    comp_frames[[i]] <- as.data.frame(as_matrix_rows(
      rows,
      8 + n_size_comp_cols[i],
      paste0("size_comp_", i)
    ))
    names(comp_frames[[i]]) <- c(
      "year", "season", "fleet", "sex", "type", "shell", "maturity",
      "nsamp", paste0("bin_", seq_len(n_size_comp_cols[i]))
    )
    pos <- pos + n_size_comp_rows[i]
  }

  growth_section <- section_lines(lines, "##_Growth data", "##_Environmental data")
  growth_nums <- numeric_lines(growth_section)
  growth_type <- as.integer(growth_nums[[1]][1])
  n_growth_obs <- as.integer(growth_nums[[2]][1])
  growth_data <- if (n_growth_obs > 0L) {
    as.data.frame(do.call(rbind, growth_nums[3:(2 + n_growth_obs)]))
  } else {
    data.frame()
  }

  env_section <- section_lines(lines, "##_Environmental data", "##_End of data file")
  env_nums <- numeric_lines(env_section)
  n_env <- as.integer(env_nums[[1]][1])
  env_ranges <- data.frame()
  env_data <- data.frame()
  if (n_env > 0L) {
    env_ranges <- as.data.frame(as_matrix_rows(env_nums[2:(1 + n_env)], 2, "env_ranges"))
    names(env_ranges) <- c("start_year", "end_year")
    n_env_data <- sum(env_ranges$end_year - env_ranges$start_year + 1L)
    env_data <- as.data.frame(as_matrix_rows(
      env_nums[(2 + n_env):(1 + n_env + n_env_data)],
      3,
      "env_data"
    ))
    names(env_data) <- c("index", "year", "value")
  }

  list(
    dimensions = setNames(as.integer(dims), names(dims)),
    n_size_sex = as.integer(n_size_sex),
    size_breaks = size_breaks,
    mid_points = size_breaks[-length(size_breaks)] + diff(size_breaks) / 2,
    season_recruitment = season_recruitment,
    season_growth = season_growth,
    season_ssb = season_ssb,
    season_N = season_N,
    m_prop_type = m_prop_type,
    m_prop = m_prop,
    fleet_names = c(fleet_names, survey_names),
    survey_names = survey_names,
    season_type = as.integer(season_type),
    catch = list(fmt = catch_fmt, rows = n_catch_rows, frames = catch_frames),
    survey = list(fmt = survey_fmt, index_type = index_type, data = survey_data),
    size_comp = list(
      fmt = comp_fmt,
      rows = n_size_comp_rows,
      cols = n_size_comp_cols,
      frames = comp_frames
    ),
    growth = list(type = growth_type, n = n_growth_obs, data = growth_data),
    environment = list(n = n_env, ranges = env_ranges, data = env_data)
  )
}

read_gmacs_pin <- function(path) {
  lines <- readLines(path, warn = FALSE)
  comment <- "^#\\s*([A-Za-z_][A-Za-z0-9_]*)\\[([0-9]+)\\].*"
  out <- list()
  names_out <- list()
  i <- 1L
  while (i <= length(lines)) {
    line <- lines[i]
    if (grepl(comment, line)) {
      group <- sub(comment, "\\1", line)
      index <- as.integer(sub(comment, "\\2", line))
      name <- trimws(sub("^.*--\\s*", "", sub(":\\s*$", "", line)))
      j <- i + 1L
      while (j <= length(lines) && !length(numeric_values(lines[j]))) {
        j <- j + 1L
      }
      value <- numeric_values(lines[j])[1]
      if (is.null(out[[group]])) {
        out[[group]] <- numeric()
        names_out[[group]] <- character()
      }
      out[[group]][index] <- value
      names_out[[group]][index] <- name
      i <- j
    }
    i <- i + 1L
  }

  collect_named_values <- function(pattern, group) {
    hits <- grep(pattern, lines)
    if (!length(hits)) {
      return()
    }
    values <- numeric(length(hits))
    labels <- character(length(hits))
    for (ii in seq_along(hits)) {
      j <- hits[ii] + 1L
      while (j <= length(lines) && !length(numeric_values(lines[j]))) {
        j <- j + 1L
      }
      values[ii] <- numeric_values(lines[j])[1]
      labels[ii] <- trimws(sub("^#\\s*", "", sub(":\\s*$", "", lines[hits[ii]])))
    }
    out[[group]] <<- values
    names_out[[group]] <<- labels
  }

  collect_single_vector <- function(label, group) {
    hit <- grep(paste0("^#\\s*", label, ":\\s*$"), lines)
    if (!length(hit)) {
      return()
    }
    j <- hit[1] + 1L
    while (j <= length(lines) && !length(numeric_values(lines[j]))) {
      j <- j + 1L
    }
    out[[group]] <<- numeric_values(lines[j])
    names_out[[group]] <<- paste0(group, "_", seq_along(out[[group]]))
  }

  collect_named_values("^#\\s*Log_fdev_", "log_fdev")
  collect_named_values("^#\\s*Log_fdov_", "log_fdov")
  collect_named_values("^#\\s*Envpar_Slx_", "slx_env_pars")
  collect_named_values("^#\\s*Slx_Devs_", "sel_devs")
  collect_single_vector("rec_ini", "rec_ini")
  collect_single_vector("rec_dev_est", "rec_dev_est")
  collect_single_vector("logit_rec_prop_est", "logit_rec_prop_est")

  for (nm in names(out)) {
    attr(out[[nm]], "parameter_names") <- names_out[[nm]]
  }
  out
}

read_bbrkc_allometry <- function(path, data) {
  lines <- readLines(path, warn = FALSE)
  dims <- data$dimensions
  nsex <- dims[["nsex"]]
  nmature <- dims[["nmature"]]
  nclass <- dims[["nclass"]]

  allometry <- numeric_lines(section_lines(lines, "##_Allometry", "##_Fecundity"))
  lw_type <- as.integer(allometry[[1]][1])
  if (lw_type != 2L) {
    stop("Only BBRKC length-weight vector inputs are currently parsed.", call. = FALSE)
  }
  mean_wt <- array(0, dim = c(nsex, nmature, nclass))
  pos <- 2L
  for (sex in seq_len(nsex)) {
    for (maturity in seq_len(nmature)) {
      mean_wt[sex, maturity, ] <- allometry[[pos]][seq_len(nclass)]
      pos <- pos + 1L
    }
  }

  fecundity <- numeric_lines(section_lines(lines, "##_Fecundity", "##_Growth parameter controls"))
  maturity <- as_matrix_rows(fecundity[seq_len(nsex)], nclass, "maturity")
  legal <- as_matrix_rows(fecundity[(nsex + 1L):(2L * nsex)], nclass, "legal")

  list(
    lw_type = lw_type,
    mean_wt = mean_wt,
    maturity = maturity,
    legal = legal,
    use_func_mat = as.integer(fecundity[[2L * nsex + 1L]][1])
  )
}

read_admb_likelihood_reference <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start <- grep("^nloglike\\s*$", lines)
  pen <- grep("^nlogPenalty\\s*$", lines)
  prior <- grep("^priorDensity\\s*$", lines)
  if (!length(start) || !length(pen) || !length(prior)) {
    stop("Could not find ADMB likelihood reference blocks.", call. = FALSE)
  }
  list(
    nloglike = unlist(numeric_lines(lines[(start[1] + 1L):(pen[1] - 1L)])),
    nlog_penalty = unlist(numeric_lines(lines[(pen[1] + 1L):(prior[1] - 1L)])),
    prior_density = unlist(numeric_lines(lines[(prior[1] + 1L):length(lines)]))
  )
}

read_admb_selectivity_block <- function(path, block) {
  lines <- readLines(path, warn = FALSE)
  start <- grep(paste0("^", block, "\\s*$"), lines)
  if (!length(start)) {
    stop("Could not find ADMB selectivity block: ", block, call. = FALSE)
  }
  possible_next <- grep("^slx_[a-z]+\\s*$", lines)
  possible_next <- possible_next[possible_next > start[1]]
  end <- if (length(possible_next)) possible_next[1] - 1L else length(lines)
  rows <- lines[(start[1] + 1L):end]
  rows <- rows[nzchar(trimws(rows))]
  rows <- rows[grepl("^\\s*[0-9]{4}\\s+", rows)]
  rows <- rows[grepl("^\\s*[0-9]{4}\\s+(male|female)\\s+", rows)]
  parsed <- lapply(rows, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    data.frame(
      year = as.integer(parts[1]),
      sex = parts[2],
      fleet = parts[3],
      size_class = seq_len(length(parts) - 3L),
      value = as.numeric(parts[-(1:3)]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parsed)
}

read_admb_fully_selected_fleet <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start <- grep("^# Fully-selected_fishing mortality by fleet\\s*$", lines)
  if (!length(start)) {
    stop("Could not find fully-selected fishing mortality by fleet block.", call. = FALSE)
  }
  next_headers <- grep("^# ", lines)
  next_headers <- next_headers[next_headers > start[1] + 1L]
  end <- if (length(next_headers)) next_headers[1] - 1L else length(lines)
  rows <- lines[(start[1] + 2L):end]
  rows <- rows[grepl("^\\s*(male|female)\\s+", rows)]
  parsed <- lapply(rows, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    values <- as.numeric(parts[-(1:3)])
    data.frame(
      sex = parts[1],
      year = as.integer(parts[2]),
      season = as.integer(parts[3]),
      fleet = seq_along(values),
      value = values,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parsed)
}

read_admb_mortality_at_size <- function(path, header) {
  lines <- readLines(path, warn = FALSE)
  start <- grep(paste0("^# ", header, "\\s*$"), lines)
  if (!length(start)) {
    stop("Could not find mortality-at-size block: ", header, call. = FALSE)
  }
  next_headers <- grep("^# ", lines)
  next_headers <- next_headers[next_headers > start[1] + 1L]
  end <- if (length(next_headers)) next_headers[1] - 1L else length(lines)
  rows <- lines[(start[1] + 2L):end]
  rows <- rows[grepl("^\\s*(male|female)\\s+", rows)]
  parsed <- lapply(rows, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    values <- as.numeric(parts[-(1:3)])
    data.frame(
      sex = parts[1],
      year = as.integer(parts[2]),
      season = as.integer(parts[3]),
      size_class = seq_along(values),
      value = values,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parsed)
}

read_admb_natural_mortality_by_class <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start <- grep("^# Natural_mortality-by-class\\s*$", lines)
  if (!length(start)) {
    stop("Could not find natural mortality by class block.", call. = FALSE)
  }
  next_headers <- grep("^# ", lines)
  next_headers <- next_headers[next_headers > start[1] + 1L]
  end <- if (length(next_headers)) next_headers[1] - 1L else length(lines)
  rows <- lines[(start[1] + 2L):end]
  rows <- rows[grepl("^\\s*[0-9]{4}\\s+(male|female)\\s+", rows)]
  parsed <- lapply(rows, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    values <- as.numeric(parts[-(1:3)])
    data.frame(
      year = as.integer(parts[1]),
      sex = parts[2],
      maturity = as.integer(parts[3]),
      size_class = seq_along(values),
      value = values,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parsed)
}

read_admb_total_mortality_by_class <- function(path, label) {
  lines <- readLines(path, warn = FALSE)
  start <- grep(paste0("^# Total mortality by size-class \\(", label, "\\)\\s*$"), lines)
  if (!length(start)) {
    stop("Could not find total mortality by class block: ", label, call. = FALSE)
  }
  next_headers <- grep("^# ", lines)
  next_headers <- next_headers[next_headers > start[1] + 1L]
  end <- if (length(next_headers)) next_headers[1] - 1L else length(lines)
  rows <- lines[(start[1] + 2L):end]
  rows <- rows[grepl("^\\s*(male|female)\\s+", rows)]
  parsed <- lapply(rows, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    values <- as.numeric(parts[-(1:4)])
    data.frame(
      sex = parts[1],
      maturity = as.integer(parts[2]),
      year = as.integer(parts[3]),
      season = as.integer(parts[4]),
      size_class = seq_along(values),
      value = values,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parsed)
}

read_admb_molt_probability <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start <- grep("^# Molt probability\\s*$", lines)
  if (!length(start)) {
    stop("Could not find molt probability block.", call. = FALSE)
  }
  next_headers <- grep("^# Growth_transition_matrix\\s*$", lines)
  end <- next_headers[next_headers > start[1]][1] - 1L
  rows <- lines[(start[1] + 2L):end]
  rows <- rows[grepl("^\\s*(male|female)\\s+", rows)]
  parsed <- lapply(rows, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    values <- as.numeric(parts[-(1:2)])
    data.frame(
      sex = parts[1],
      year = as.integer(parts[2]),
      size_class = seq_along(values),
      value = values,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parsed)
}

read_admb_growth_transition <- function(path, nclass = 20L) {
  lines <- readLines(path, warn = FALSE)
  headers <- grep("^#growth_matrix for \\(sex, increment_no\\):", lines)
  if (!length(headers)) {
    stop("Could not find growth transition matrix blocks.", call. = FALSE)
  }
  parsed <- vector("list", length(headers))
  for (i in seq_along(headers)) {
    header <- strsplit(trimws(lines[headers[i]]), "\\s+")[[1]]
    sex <- header[length(header) - 1L]
    increment_no <- as.integer(header[length(header)])
    matrix_lines <- lines[(headers[i] + 1L):(headers[i] + nclass)]
    mat <- as_matrix_rows(numeric_lines(matrix_lines), nclass, "growth_transition")
    parsed[[i]] <- data.frame(
      sex = sex,
      increment_no = increment_no,
      to_size_class = rep(seq_len(nclass), each = nclass),
      from_size_class = rep(seq_len(nclass), times = nclass),
      value = as.vector(t(mat)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, parsed)
}

read_admb_overall_summary <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start <- grep("^#Overall_summary\\s*$", lines)
  if (!length(start)) {
    stop("Could not find overall summary block.", call. = FALSE)
  }
  header <- strsplit(trimws(sub("^#", "", lines[start[1] + 1L])), "\\s+")[[1]]
  next_headers <- grep("^#[-A-Za-z_]", lines)
  next_headers <- next_headers[next_headers > start[1] + 1L]
  end <- if (length(next_headers)) next_headers[1] - 1L else length(lines)
  rows <- lines[(start[1] + 2L):end]
  rows <- rows[grepl("^\\s*[0-9]{4}\\s+", rows)]
  values <- as_matrix_rows(numeric_lines(rows), length(header), "overall_summary")
  out <- as.data.frame(values)
  names(out) <- make.names(header, unique = TRUE)
  names(out)[1] <- "Year"
  out
}

read_admb_numbers_at_size <- function(path, block) {
  lines <- readLines(path, warn = FALSE)
  start <- grep(paste0("^# ", gsub("([()])", "\\\\\\1", block), "\\s*$"), lines)
  if (!length(start)) {
    stop("Could not find ADMB numbers-at-size block: ", block, call. = FALSE)
  }
  next_headers <- grep("^# N\\(", lines)
  next_headers <- next_headers[next_headers > start[1]]
  end <- if (length(next_headers)) next_headers[1] - 1L else length(lines)
  rows <- lines[(start[1] + 2L):end]
  rows <- rows[grepl("^\\s*[0-9]{4}\\s+", rows)]
  parsed <- lapply(rows, function(x) {
    values <- numeric_values(x)
    data.frame(
      year = as.integer(values[1]),
      size_class = seq_along(values[-1]),
      value = values[-1],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parsed)
}

read_admb_catch_fit_summary <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start <- grep("^#Catch_fit_summary\\s*$", lines)
  if (!length(start)) {
    stop("Could not find catch fit summary block.", call. = FALSE)
  }
  next_headers <- grep("^#--------------------------------------------------------------------------------------------", lines)
  next_headers <- next_headers[next_headers > start[1]]
  end <- if (length(next_headers)) next_headers[1] - 1L else length(lines)
  rows <- lines[(start[1] + 2L):end]
  rows <- rows[grepl("^\\s*[0-9]+\\s+", rows)]
  parsed <- lapply(rows, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    data.frame(
      series = as.integer(parts[1]),
      year = as.integer(parts[2]),
      fleet = parts[3],
      season = as.integer(parts[4]),
      sex = tolower(parts[5]),
      obs = as.numeric(parts[6]),
      cv = as.numeric(parts[7]),
      type = tolower(parts[8]),
      units = tolower(parts[9]),
      mult = as.numeric(parts[10]),
      effort = as.numeric(parts[11]),
      discard_mortality = as.numeric(parts[12]),
      predicted = as.numeric(parts[13]),
      residual = as.numeric(parts[14]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parsed)
}

read_admb_size_fit_summary <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start <- grep("^Size_fit_summary\\s*$", lines)
  if (!length(start)) {
    stop("Could not find size fit summary block.", call. = FALSE)
  }
  next_headers <- grep("^#Size data: standard deviation and median", lines)
  next_headers <- next_headers[next_headers > start[1]]
  end <- if (length(next_headers)) next_headers[1] - 1L else length(lines)
  rows <- lines[(start[1] + 2L):end]
  rows <- rows[grepl("^\\s*[0-9]+\\s+", rows)]
  parsed <- vector("list", length(rows))
  for (i in seq_along(rows)) {
    parts <- strsplit(trimws(rows[i]), "\\s+")[[1]]
    values <- as.numeric(parts[-(1:10)])
    nbin <- length(values) / 2L
    parsed[[i]] <- rbind(
      data.frame(
        original_series = as.integer(parts[1]),
        modified_series = as.integer(parts[2]),
        year = as.integer(parts[3]),
        fleet = parts[4],
        season = as.integer(parts[5]),
        sex = tolower(parts[6]),
        type = tolower(parts[7]),
        shell = tolower(parts[8]),
        maturity = tolower(parts[9]),
        nsamp = as.numeric(parts[10]),
        vector = "obs",
        size_bin = seq_len(nbin),
        value = values[seq_len(nbin)],
        stringsAsFactors = FALSE
      ),
      data.frame(
        original_series = as.integer(parts[1]),
        modified_series = as.integer(parts[2]),
        year = as.integer(parts[3]),
        fleet = parts[4],
        season = as.integer(parts[5]),
        sex = tolower(parts[6]),
        type = tolower(parts[7]),
        shell = tolower(parts[8]),
        maturity = tolower(parts[9]),
        nsamp = as.numeric(parts[10]),
        vector = "pred",
        size_bin = seq_len(nbin),
        value = values[nbin + seq_len(nbin)],
        stringsAsFactors = FALSE
      )
    )
  }
  do.call(rbind, parsed)
}

read_admb_index_fit_summary <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start <- grep("^#Index_fit_summary\\s*$", lines)
  if (!length(start)) {
    stop("Could not find index fit summary block.", call. = FALSE)
  }
  next_headers <- grep("^# CPUE: standard deviation and median", lines)
  next_headers <- next_headers[next_headers > start[1]]
  end <- if (length(next_headers)) next_headers[1] - 1L else length(lines)
  rows <- lines[(start[1] + 2L):end]
  rows <- rows[grepl("^\\s*[0-9]+\\s+", rows)]
  parsed <- lapply(rows, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    data.frame(
      series = as.integer(parts[1]),
      year = as.integer(parts[2]),
      fleet = parts[3],
      season = as.integer(parts[4]),
      sex = tolower(parts[5]),
      maturity = tolower(parts[6]),
      obs = as.numeric(parts[7]),
      base_cv = as.numeric(parts[8]),
      actual_cv = as.numeric(parts[9]),
      units = tolower(parts[10]),
      q = as.numeric(parts[11]),
      time = as.numeric(parts[12]),
      predicted = as.numeric(parts[13]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parsed)
}

read_bbrkc_control_phases <- function(path) {
  lines <- readLines(path, warn = FALSE)

  numeric_rows <- function(start_pattern, end_pattern, ncol) {
    rows <- numeric_lines(section_lines(lines, start_pattern, end_pattern))
    rows[lengths(rows) == ncol]
  }

  theta <- numeric_rows("##_Key parameter controls", "##_Allometry", 7)
  growth <- numeric_rows("#_Growth increment model controls", "#_Molt probability controls", 7)
  molt <- numeric_rows("#_Molt probability controls", "#_Mature probability controls", 7)
  mature <- numeric_rows("#_Mature probability controls", "#_Custom growth-increment", 7)
  selectivity <- numeric_rows("#_Selectivity parameter controls", "#_Retention parameter controls", 19)
  retention <- numeric_rows("#_Retention parameter controls", "#_Number of asymptotic retention parameter", 19)
  asymret <- numeric_rows("#_Asymptotic parameter controls", "#_Environmental parameters Control", 7)
  fishing <- numeric_rows("## Penalties for the average fishing mortality rate", "## Size composition data control", 12)
  size_comp <- numeric_lines(section_lines(lines, "## Size composition data control", "## Time-varying Natural mortality controls"))
  mdev <- numeric_rows("## Time-varying Natural mortality controls", "## Tagging controls", 5)
  mmat <- numeric_rows("##  Immature/mature natural mortality", "## Other (additional) controls", 7)
  q <- numeric_rows("## Priors for catchability", "## Additional CV controls", 10)
  add_cv <- numeric_rows("## Additional CV controls", "## Penalties for the average fishing mortality rate", 7)
  add_cv_section <- numeric_lines(section_lines(lines, "## Additional CV controls", "## Penalties for the average fishing mortality rate"))

  growth_section <- numeric_lines(section_lines(lines, "##_Growth parameter controls", "#_Selectivity"))
  b_use_custom_growth_matrix <- as.integer(growth_section[[1]][1])
  b_use_growth_increment_model <- as.integer(growth_section[[2]][1])
  b_use_custom_molt_probability <- as.integer(growth_section[[3]][1])
  max_recruit_size <- as.integer(growth_section[[4]])
  n_size_inc_varies <- as.integer(growth_section[[5]])
  size_change_rows <- if (any(n_size_inc_varies > 1L)) sum(n_size_inc_varies > 1L) else 0L
  size_change_start <- 6L
  size_change_end <- size_change_start + size_change_rows - 1L
  iyrs_size_changes <- vector("list", length(n_size_inc_varies))
  if (size_change_rows > 0L) {
    change_values <- growth_section[size_change_start:size_change_end]
    pos <- 1L
    for (sex in seq_along(n_size_inc_varies)) {
      n_change <- n_size_inc_varies[sex] - 1L
      if (n_change > 0L) {
        iyrs_size_changes[[sex]] <- as.integer(change_values[[pos]])
        pos <- pos + 1L
      } else {
        iyrs_size_changes[[sex]] <- integer()
      }
    }
  } else {
    iyrs_size_changes <- lapply(n_size_inc_varies, function(x) integer())
  }
  n_molt_varies_pos <- size_change_end + 1L
  n_molt_varies <- as.integer(growth_section[[n_molt_varies_pos]])
  molt_change_rows <- if (any(n_molt_varies > 1L)) sum(n_molt_varies > 1L) else 0L
  molt_change_start <- n_molt_varies_pos + 1L
  molt_change_end <- molt_change_start + molt_change_rows - 1L
  iyrs_molt_changes <- vector("list", length(n_molt_varies))
  if (molt_change_rows > 0L) {
    change_values <- growth_section[molt_change_start:molt_change_end]
    pos <- 1L
    for (sex in seq_along(n_molt_varies)) {
      n_change <- n_molt_varies[sex] - 1L
      if (n_change > 0L) {
        iyrs_molt_changes[[sex]] <- as.integer(change_values[[pos]])
        pos <- pos + 1L
      } else {
        iyrs_molt_changes[[sex]] <- integer()
      }
    }
  } else {
    iyrs_molt_changes <- lapply(n_molt_varies, function(x) integer())
  }
  beta_par_relative <- as.integer(growth_section[[molt_change_end + 1L]][1])

  other <- numeric_lines(section_lines(lines, "## Other (additional) controls", "## Emphasis factor"))
  mdev_section <- numeric_lines(section_lines(lines, "## Time-varying Natural mortality controls", "## Tagging controls"))
  m_type <- as.integer(mdev_section[[1]][1])
  m_rel_female <- as.integer(mdev_section[[2]][1])
  mdev_phase <- mdev_section[[3]][1]
  mdev_sd <- mdev_section[[4]][1]
  m_n_nodes_sex <- as.integer(unlist(mdev_section[5:6], use.names = FALSE))
  m_nodeyear_sex <- as_matrix_rows(
    mdev_section[7:(6 + length(m_n_nodes_sex))],
    max(m_n_nodes_sex),
    "m_nodeyear_sex"
  )
  n_size_devs_pos <- 7 + length(m_n_nodes_sex)
  n_size_devs <- as.integer(mdev_section[[n_size_devs_pos]][1])
  m_size_nodeyear <- if (n_size_devs > 0L) {
    as.integer(mdev_section[[n_size_devs_pos + 1L]])
  } else {
    integer()
  }
  mdev_specific_start <- n_size_devs_pos + 1L + as.integer(n_size_devs > 0L)
  mdev_specific_ival <- as.integer(mdev_section[[mdev_specific_start]][1])
  mdev_control_start <- mdev_specific_start + 1L
  mdev_controls <- as_matrix_rows(
    mdev_section[mdev_control_start:(mdev_control_start + length(mdev) - 1L)],
    5,
    "mdev_controls"
  )

  n_mdev_par_count <- switch(
    as.character(m_type),
    "0" = rep(0L, length(m_n_nodes_sex)),
    "1" = rep(as.integer(numeric_values(line_after(lines, "#_End year")) -
      numeric_values(line_after(lines, "#_Start year"))), length(m_n_nodes_sex)),
    "2" = m_n_nodes_sex,
    "3" = m_n_nodes_sex,
    "4" = m_n_nodes_sex / 2L,
    "5" = m_n_nodes_sex,
    "6" = m_n_nodes_sex,
    stop("Unsupported natural mortality type: ", m_type, call. = FALSE)
  )
  if (length(n_mdev_par_count) > 1L && m_type == 1L) {
    n_mdev_par_count[2] <- n_mdev_par_count[1]
  }

  mmat_section <- numeric_lines(section_lines(lines, "##  Immature/mature natural mortality", "## Other (additional) controls"))
  m_maturity <- as.integer(mmat_section[[1]][1])
  rec_control <- other[seq_len(12)]
  init_sex_ratio <- rec_control[[6]][1]
  emphasis <- numeric_lines(section_lines(lines, "## Emphasis factor", "## End of control file"))
  catch_emphasis <- emphasis[[1]]
  fdev_penalty <- as_matrix_rows(emphasis[2:(1 + length(fishing))], 4, "fdev_penalty")
  penalty_emphasis <- unlist(
    emphasis[(2 + length(fishing)):(14 + length(fishing))],
    use.names = FALSE
  )

  log_vn_phase_in <- size_comp[[4]]
  comp_aggregator <- size_comp[[5]]
  log_vn_phase <- rep(NA_real_, max(comp_aggregator))
  for (i in seq_along(log_vn_phase_in)) {
    log_vn_phase[comp_aggregator[i]] <- log_vn_phase_in[i]
  }

  list(
    theta = vapply(theta, `[`, numeric(1), 4),
    Grwth = c(vapply(growth, `[`, numeric(1), 4), vapply(molt, `[`, numeric(1), 4), vapply(mature, `[`, numeric(1), 4)),
    log_slx_pars = c(vapply(selectivity, `[`, numeric(1), 11), vapply(retention, `[`, numeric(1), 11)),
    Asymret = vapply(asymret, `[`, numeric(1), 7),
    slx_env_pars = -1,
    sel_devs = size_comp[[1]][1] * 0 - 1,
    log_fbar = vapply(fishing, `[`, numeric(1), 5),
    log_foff = vapply(fishing, `[`, numeric(1), 6),
    log_vn = log_vn_phase,
    survey_q = vapply(q, `[`, numeric(1), 4),
    log_add_cv = vapply(add_cv, `[`, numeric(1), 4),
    m_dev_est = vapply(mdev, `[`, numeric(1), 4),
    m_mat_mult = vapply(mmat, `[`, numeric(1), 4),
    growth_controls = list(
      b_use_custom_growth_matrix = b_use_custom_growth_matrix,
      b_use_growth_increment_model = b_use_growth_increment_model,
      b_use_custom_molt_probability = b_use_custom_molt_probability,
      max_recruit_size = max_recruit_size,
      n_size_inc_varies = n_size_inc_varies,
      iyrs_size_changes = iyrs_size_changes,
      n_molt_varies = n_molt_varies,
      iyrs_molt_changes = iyrs_molt_changes,
      beta_par_relative = beta_par_relative
    ),
    size_comp_controls = list(
      likelihood_type_in = as.integer(size_comp[[1]]),
      tail_compression_in = as.integer(size_comp[[2]]),
      nvn_ival_in = size_comp[[3]],
      nvn_phz_in = as.integer(size_comp[[4]]),
      comp_aggregator = as.integer(size_comp[[5]]),
      lf_catch_in = as.integer(size_comp[[6]]),
      lambda = size_comp[[7]],
      emphasis = size_comp[[8]]
    ),
    likelihood_emphasis = list(
      catch = catch_emphasis,
      cpue_lambda = vapply(q, `[`, numeric(1), 9),
      cpue_emphasis = vapply(q, `[`, numeric(1), 10),
      fdev_penalty = fdev_penalty,
      penalty = penalty_emphasis
    ),
    prior_controls = list(
      theta = do.call(rbind, theta),
      Grwth = do.call(rbind, c(growth, molt, mature)),
      log_slx_pars = do.call(rbind, c(selectivity, retention)),
      Asymret = do.call(rbind, asymret),
      m_dev_est = do.call(rbind, mdev),
      q = do.call(rbind, q),
      add_cv = do.call(rbind, add_cv),
      m_mat_mult = do.call(rbind, mmat)
    ),
    q_controls = list(
      q_ival = vapply(q, `[`, numeric(1), 1),
      q_lb = vapply(q, `[`, numeric(1), 2),
      q_ub = vapply(q, `[`, numeric(1), 3),
      q_phz = vapply(q, `[`, numeric(1), 4),
      prior_qtype = vapply(q, `[`, numeric(1), 5),
      prior_p1 = vapply(q, `[`, numeric(1), 6),
      prior_p2 = vapply(q, `[`, numeric(1), 7),
      q_anal = as.integer(vapply(q, `[`, numeric(1), 8)),
      cpue_lambda = vapply(q, `[`, numeric(1), 9),
      cpue_emphasis = vapply(q, `[`, numeric(1), 10),
      add_cv_links = as.integer(add_cv_section[[length(add_cv_section)]])
    ),
    rec_controls = list(
      rdv_syr = as.integer(rec_control[[1]][1]),
      rdv_eyr = as.integer(rec_control[[2]][1]),
      terminal_molt = as.integer(rec_control[[3]][1]),
      rdv_phz = as.integer(rec_control[[4]][1]),
      rec_prop_phz = as.integer(rec_control[[5]][1]),
      init_sex_ratio = init_sex_ratio,
      init_logit_sex_ratio = -log((1 - init_sex_ratio) / init_sex_ratio),
      rec_ini_phz = as.integer(rec_control[[7]][1]),
      initialize_unfished = as.integer(rec_control[[8]][1]),
      spr_lambda = rec_control[[9]][1],
      stock_recruit_flag = as.integer(rec_control[[10]][1]),
      brp_rec_sex_ratio = as.integer(rec_control[[11]][1]),
      n_year_equilibrium = as.integer(rec_control[[12]][1])
    ),
    m_controls = list(
      m_type = m_type,
      m_rel_female = m_rel_female,
      mdev_phase = mdev_phase,
      mdev_sd = mdev_sd,
      m_n_nodes_sex = m_n_nodes_sex,
      m_nodeyear_sex = m_nodeyear_sex,
      n_size_devs = n_size_devs,
      m_size_nodeyear = m_size_nodeyear,
      mdev_specific_ival = mdev_specific_ival,
      mdev_spec = as.integer(mdev_controls[, 5]),
      n_mdev_par_count = as.integer(n_mdev_par_count),
      m_maturity = m_maturity
    ),
    rec_ini = rep(other[[7]][1], length.out = 0),
    rdv_phz = other[[4]][1],
    rec_prop_phz = other[[5]][1],
    rec_ini_phz = other[[7]][1]
  )
}

read_bbrkc_inputs <- function(root = "build/BBRKC") {
  main <- read_gmacs_main(file.path(root, "gmacs.dat"))
  data <- read_gmacs_data(file.path(root, main$data_file))
  control_path <- file.path(root, main$control_file)
  data$allometry <- read_bbrkc_allometry(control_path, data)
  list(
    main = main,
    data = data,
    control_phases = read_bbrkc_control_phases(control_path),
    parameters = read_gmacs_pin(file.path(root, "gmacs.pin")),
    admb_reference = read_admb_likelihood_reference(file.path(root, "Gmacsall.out")),
    admb_selectivity = list(
      capture = read_admb_selectivity_block(file.path(root, "gmacs.rep"), "slx_capture"),
      retained = read_admb_selectivity_block(file.path(root, "gmacs.rep"), "slx_retaind"),
      discard = read_admb_selectivity_block(file.path(root, "gmacs.rep"), "slx_discard")
    ),
    admb_fishing_mortality = list(
      fully_selected_fleet = read_admb_fully_selected_fleet(file.path(root, "Gmacsall.out")),
      at_size_continuous = read_admb_mortality_at_size(
        file.path(root, "Gmacsall.out"),
        "Fishing mortality-at-size by sex and season \\(Continuous\\)"
      )
    ),
    admb_natural_mortality = read_admb_natural_mortality_by_class(file.path(root, "Gmacsall.out")),
    admb_total_mortality = list(
      continuous = read_admb_total_mortality_by_class(file.path(root, "Gmacsall.out"), "continuous"),
      discrete = read_admb_total_mortality_by_class(file.path(root, "Gmacsall.out"), "discrete")
    ),
    admb_growth = list(
      molt_probability = read_admb_molt_probability(file.path(root, "Gmacsall.out")),
      growth_transition = read_admb_growth_transition(
        file.path(root, "Gmacsall.out"),
        nclass = data$dimensions[["nclass"]]
      )
    ),
    admb_numbers_at_size = list(
      total = read_admb_numbers_at_size(file.path(root, "Gmacsall.out"), "N(total)"),
      males = read_admb_numbers_at_size(file.path(root, "Gmacsall.out"), "N(males)"),
      females = read_admb_numbers_at_size(file.path(root, "Gmacsall.out"), "N(females)"),
      males_new = read_admb_numbers_at_size(file.path(root, "Gmacsall.out"), "N(males_new)"),
      females_new = read_admb_numbers_at_size(file.path(root, "Gmacsall.out"), "N(females_new)"),
      males_old = read_admb_numbers_at_size(file.path(root, "Gmacsall.out"), "N(males_old)"),
      females_old = read_admb_numbers_at_size(file.path(root, "Gmacsall.out"), "N(females_old)")
    ),
    admb_catch_fit = read_admb_catch_fit_summary(file.path(root, "Gmacsall.out")),
    admb_index_fit = read_admb_index_fit_summary(file.path(root, "Gmacsall.out")),
    admb_size_fit = read_admb_size_fit_summary(file.path(root, "Gmacsall.out")),
    admb_overall_summary = read_admb_overall_summary(file.path(root, "Gmacsall.out")),
    paths = list(
      root = root,
      data = file.path(root, main$data_file),
      control = file.path(root, main$control_file),
      projection = file.path(root, main$project_file),
      pin = file.path(root, "gmacs.pin"),
      admb_output = file.path(root, "Gmacsall.out")
    )
  )
}

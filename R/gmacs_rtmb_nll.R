dnorm_nll <- function(x, mean = 0, sd = 1) {
  -sum(dnorm(x, mean = mean, sd = sd, log = TRUE))
}

admb_dnorm_nll <- function(x, sd, mean = 0) {
  sum(0.5 * log(2 * pi) + log(sd) + 0.5 * ((x - mean) / sd)^2)
}

admb_prior_nll <- function(type, x, p1, p2) {
  switch(
    as.character(as.integer(type)),
    "0" = log(p2 - p1),
    "1" = admb_dnorm_nll(x, p2, p1),
    "2" = 0.5 * log(2 * pi) + log(p2) + log(x) + 0.5 * ((log(x) - log(p1)) / p2)^2,
    "3" = -(lgamma(p1 + p2) - lgamma(p1) - lgamma(p2) +
      (p1 - 1) * log(x) + (p2 - 1) * log(1 - x)),
    "4" = -(p1 * log(p2) - lgamma(p1) + (p1 - 1) * log(x) - p2 * x),
    stop("Unsupported prior type: ", type, call. = FALSE)
  )
}

first_difference_nll <- function(x, sd = 1) {
  if (length(x) < 2L) {
    return(sum(x * 0))
  }
  out <- 0
  for (i in 2:length(x)) {
    out <- out + admb_dnorm_nll(x[i] - x[i - 1L], sd)
  }
  out
}

logistic_selectivity <- function(size, x50, slope) {
  1 / (1 + exp(-(size - x50) / slope))
}

logistic95_selectivity <- function(size, x50, x95) {
  1 / (1 + exp(-log(19) * ((size - x50) / (x95 - x50))))
}

plogis_gmacs <- function(x, location, scale) {
  1 / (1 + exp(-(x - location) / scale))
}

calc_bbrkc_capture_selectivity <- function(log_slx_pars, data) {
  years <- data$years
  n_years <- length(years)
  n_class <- length(data$mid_points)
  out <- array(0, dim = c(6, 2, n_years, n_class))

  fill_logistic <- function(fleet, sex, year_min, year_max, p_index) {
    rows <- years >= year_min & years <= year_max
    sel <- logistic_selectivity(
      data$mid_points,
      exp(log_slx_pars[p_index]),
      exp(log_slx_pars[p_index + 1L])
    )
    sel <- sel / sel[length(sel)]
    log_sel <- log(sel)
    for (iy in which(rows)) {
      out[fleet, sex, iy, ] <<- log_sel
    }
  }

  fill_logistic(1, 1, 1975, 2020, 1)
  fill_logistic(1, 2, 1975, 2020, 3)
  fill_logistic(2, 1, 1975, 2020, 5)
  fill_logistic(2, 2, 1975, 2020, 5)
  fill_logistic(3, 1, 1975, 2020, 7)
  fill_logistic(3, 2, 1975, 2020, 9)
  fill_logistic(4, 1, 1975, 2020, 11)
  fill_logistic(4, 2, 1975, 2020, 11)
  fill_logistic(5, 1, 1975, 1981, 13)
  fill_logistic(5, 2, 1975, 1981, 13)
  fill_logistic(5, 1, 1982, 2021, 15)
  fill_logistic(5, 2, 1982, 2021, 15)
  fill_logistic(6, 1, 1975, 2021, 17)
  fill_logistic(6, 2, 1975, 2021, 17)

  ## BBRKC control file sets NMFS_Trawl selectivity as embedded within BSFRF.
  out[5, , , ] <- out[5, , , ] + out[6, , , ]

  out
}

calc_bbrkc_retention_selectivity <- function(log_slx_pars, Asymret, data) {
  years <- data$years
  n_years <- length(years)
  n_class <- length(data$mid_points)
  log_retained <- array(log(0), dim = c(6, 2, n_years, n_class))
  log_discard <- array(log(1), dim = c(6, 2, n_years, n_class))

  fill_pot_male <- function(year_min, year_max, p_index) {
    rows <- years >= year_min & years <= year_max
    ret <- logistic_selectivity(
      data$mid_points,
      exp(log_slx_pars[p_index]),
      exp(log_slx_pars[p_index + 1L])
    )
    ret <- ret / ret[length(ret)]
    for (iy in which(rows)) {
      year <- years[iy]
      high_grade <- if (year == 1975) 1 - Asymret[1] else 1
      retained <- ret * high_grade
      log_retained[1, 1, iy, ] <<- log(retained)
      log_discard[1, 1, iy, ] <<- log(1 - retained + 1e-08)
    }
  }

  fill_pot_male(1975, 2004, 19)
  fill_pot_male(2005, 2020, 21)

  list(
    log_slx_retaind = log_retained,
    log_slx_discard = log_discard
  )
}

split_parameter_by_fleet <- function(values, indices) {
  lapply(indices, function(i) values[i])
}

build_bbrkc_fishing_hits <- function(data) {
  dims <- data$dimensions
  years <- data$years
  n_years <- length(years)
  nseason <- dims[["nseason"]]
  nfleet <- dims[["nfleet"]]

  fhit <- array(0, dim = c(n_years, nseason, nfleet))
  yhit <- array(0, dim = c(n_years, nseason, nfleet))
  dmr <- matrix(0, nrow = n_years, ncol = nfleet)

  for (frame in data$catch$frames) {
    for (i in seq_len(nrow(frame))) {
      year <- as.integer(frame$year[i])
      if (year < dims[["syr"]] || year > dims[["nyr"]]) {
        next
      }
      iy <- match(year, years)
      season <- as.integer(frame$season[i])
      fleet <- as.integer(frame$fleet[i])
      sex <- as.integer(frame$sex[i])
      if (fhit[iy, season, fleet] == 0) {
        fhit[iy, season, fleet] <- 1
        dmr[iy, fleet] <- frame$discard_mortality[i]
      }
      if (yhit[iy, season, fleet] == 0 && sex == 2L) {
        yhit[iy, season, fleet] <- 1
        dmr[iy, fleet] <- frame$discard_mortality[i]
      }
    }
  }

  list(fhit = fhit, yhit = yhit, dmr = dmr)
}

calc_bbrkc_fishing_mortality <- function(par, data, model) {
  dims <- data$dimensions
  years <- data$years
  n_years <- length(years)
  nfleet <- dims[["nfleet"]]
  nsex <- dims[["nsex"]]
  nseason <- dims[["nseason"]]
  nclass <- dims[["nclass"]]

  hits <- build_bbrkc_fishing_hits(data)
  log_fdev <- split_parameter_by_fleet(par$log_fdev, data$fdev_indices)
  log_fdov <- split_parameter_by_fleet(par$log_fdov, data$fdov_indices)

  ft <- array(0, dim = c(nfleet, nsex, n_years, nseason))
  fout <- matrix(0, nrow = nfleet, ncol = n_years)
  F <- array(0, dim = c(nsex, n_years, nseason, nclass))
  F2 <- array(1e-10, dim = c(nsex, n_years, nseason, nclass))
  slx_nret <- matrix(0, nrow = nsex, ncol = nfleet)
  slx_nret[1, 1] <- 1

  for (fleet in seq_len(nfleet)) {
    for (sex in seq_len(nsex)) {
      ik <- 1L
      yk <- 1L
      for (iy in seq_len(n_years)) {
        for (season in seq_len(nseason)) {
          if (hits$fhit[iy, season, fleet] > 0) {
            log_ftmp <- par$log_fbar[fleet] + log_fdev[[fleet]][ik]
            ik <- ik + 1L
            fout[fleet, iy] <- exp(log_ftmp)
            if (sex == 2L) {
              log_ftmp <- log_ftmp + par$log_foff[fleet]
            }
            if (sex == 2L && hits$yhit[iy, season, fleet] > 0) {
              log_ftmp <- log_ftmp + log_fdov[[fleet]][yk]
              yk <- yk + 1L
            }
            ft[fleet, sex, iy, season] <- exp(log_ftmp)
            xi <- hits$dmr[iy, fleet]
            sel <- exp(model$log_slx_capture[fleet, sex, iy, ]) + 1e-10
            ret <- exp(model$log_slx_retaind[fleet, sex, iy, ]) * slx_nret[sex, fleet]
            vul <- sel * (ret + (1 - ret) * xi)
            F[sex, iy, season, ] <- F[sex, iy, season, ] + ft[fleet, sex, iy, season] * vul
            F2[sex, iy, season, ] <- F2[sex, iy, season, ] + ft[fleet, sex, iy, season] * sel
          }
        }
      }
    }
  }

  list(
    fhit = hits$fhit,
    yhit = hits$yhit,
    dmr = hits$dmr,
    ft = ft,
    fout = fout,
    F = F,
    F2 = F2
  )
}

calc_bbrkc_natural_mortality <- function(par, data, model) {
  dims <- data$dimensions
  controls <- data$m_controls
  years <- data$years
  n_years <- length(years)
  nsex <- dims[["nsex"]]
  nmature <- dims[["nmature"]]
  nclass <- dims[["nclass"]]

  m_dev <- par$m_dev_est
  if (controls$n_size_devs > 0L) {
    m_dev <- par$m_dev_est[seq_len(length(par$m_dev_est) - controls$n_size_devs)]
  }
  for (i in seq_along(m_dev)) {
    spec <- controls$mdev_spec[i]
    if (spec < 0L) {
      m_dev[i] <- par$m_dev_est[-spec]
    }
  }

  m_dev_sex <- vector("list", nsex)
  pos <- 1L
  for (sex in seq_len(nsex)) {
    n_dev <- controls$n_mdev_par_count[sex]
    m_dev_sex[[sex]] <- m_dev[pos:(pos + n_dev - 1L)]
    pos <- pos + n_dev
  }

  m_mult <- rep(1, nclass)
  if (controls$n_size_devs > 0L) {
    size_dev_offset <- length(par$m_dev_est) - controls$n_size_devs
    for (i in seq_len(controls$n_size_devs)) {
      m_mult[controls$m_size_nodeyear[i]:nclass] <- par$m_dev_est[size_dev_offset + i]
    }
  }

  M <- array(0, dim = c(nsex, nmature, n_years, nclass))
  for (maturity in seq_len(nmature)) {
    for (sex in seq_len(nsex)) {
      annual <- rep(model$M0[sex], n_years)
      if (controls$m_maturity == 1L && maturity == 2L) {
        annual <- rep(model$M0[sex] * exp(par$m_mat_mult[sex]), n_years)
      }

      if (controls$m_type == 6L) {
        for (idev in seq_len(controls$n_mdev_par_count[sex] - 1L)) {
          from <- controls$m_nodeyear_sex[sex, idev]
          to <- controls$m_nodeyear_sex[sex, idev + 1L] - 1L
          rows <- years >= from & years <= to
          annual[rows] <- annual[rows] * exp(m_dev_sex[[sex]][idev])
        }
      } else if (controls$m_type != 0L) {
        stop("Natural mortality type ", controls$m_type, " is not ported yet.", call. = FALSE)
      }

      for (iy in seq_len(n_years)) {
        M[sex, maturity, iy, ] <- annual[iy] * m_mult
      }
    }
  }

  M
}

calc_bbrkc_total_mortality <- function(data, model) {
  dims <- data$dimensions
  nsex <- dims[["nsex"]]
  nmature <- dims[["nmature"]]
  nseason <- dims[["nseason"]]
  nclass <- dims[["nclass"]]
  n_years <- length(data$years)

  Z <- array(0, dim = c(nsex, nmature, n_years, nseason, nclass))
  Z2 <- array(0, dim = c(nsex, nmature, n_years, nseason, nclass))
  S <- array(0, dim = c(nsex, nmature, n_years, nseason, nclass))

  for (sex in seq_len(nsex)) {
    for (maturity in seq_len(nmature)) {
      for (iy in seq_len(n_years)) {
        for (season in seq_len(nseason)) {
          Z[sex, maturity, iy, season, ] <-
            data$m_prop[iy, season] * model$M[sex, maturity, iy, ] +
            model$F[sex, iy, season, ]
          Z2[sex, maturity, iy, season, ] <-
            data$m_prop[iy, season] * model$M[sex, maturity, iy, ] +
            model$F2[sex, iy, season, ]
          if (data$season_type[season] == 0L) {
            S[sex, maturity, iy, season, ] <- 1 -
              Z[sex, maturity, iy, season, ] / Z2[sex, maturity, iy, season, ] *
                (1 - exp(-Z2[sex, maturity, iy, season, ]))
          } else {
            S[sex, maturity, iy, season, ] <- exp(-Z[sex, maturity, iy, season, ])
          }
        }
      }
    }
  }

  list(Z = Z, Z2 = Z2, S = S)
}

unpack_bbrkc_growth_parameters <- function(par, data) {
  controls <- data$growth_controls
  nsex <- data$dimensions[["nsex"]]
  nclass <- data$dimensions[["nclass"]]
  max_blocks <- max(controls$n_size_inc_varies)
  molt_increment <- array(0, dim = c(nsex, max_blocks, nclass))
  gscale <- matrix(0, nrow = nsex, ncol = max_blocks)
  molt_mu <- matrix(0, nrow = nsex, ncol = max(controls$n_molt_varies))
  molt_cv <- matrix(0, nrow = nsex, ncol = max(controls$n_molt_varies))

  pos <- 1L
  for (sex in seq_len(nsex)) {
    for (block in seq_len(controls$n_size_inc_varies[sex])) {
      molt_increment[sex, block, ] <- par$Grwth[pos:(pos + nclass - 1L)]
      pos <- pos + nclass
      if (controls$beta_par_relative == 1L && block > 1L) {
        gscale[sex, block] <- exp(par$Grwth[pos]) * gscale[sex, 1L]
      } else {
        gscale[sex, block] <- par$Grwth[pos]
      }
      pos <- pos + 1L
    }
  }
  for (sex in seq_len(nsex)) {
    for (block in seq_len(controls$n_molt_varies[sex])) {
      molt_mu[sex, block] <- par$Grwth[pos]
      molt_cv[sex, block] <- par$Grwth[pos + 1L]
      pos <- pos + 2L
    }
  }

  list(
    molt_increment = molt_increment,
    gscale = gscale,
    molt_mu = molt_mu,
    molt_cv = molt_cv
  )
}

calc_bbrkc_molting_probability <- function(growth, data) {
  controls <- data$growth_controls
  years <- data$years
  n_years <- length(years)
  nsex <- data$dimensions[["nsex"]]
  nclass <- data$dimensions[["nclass"]]
  out <- array(0, dim = c(nsex, n_years, nclass))

  if (controls$b_use_custom_molt_probability != 2L) {
    stop("Only logistic BBRKC molt probability is ported.", call. = FALSE)
  }

  for (sex in seq_len(nsex)) {
    for (block in seq_len(controls$n_molt_varies[sex])) {
      rows <- rep(TRUE, n_years)
      if (block > 1L) {
        rows <- years >= controls$iyrs_molt_changes[[sex]][block - 1L]
      }
      mu <- growth$molt_mu[sex, block]
      sd <- mu * growth$molt_cv[sex, block]
      values <- 1 - plogis_gmacs(data$mid_points, mu, sd)
      for (iy in which(rows)) {
        out[sex, iy, ] <- values
      }
    }
  }

  out
}

calc_bbrkc_growth_transition <- function(growth, data) {
  if (!is.null(data$growth_transition_fixed)) {
    return(data$growth_transition_fixed)
  }

  controls <- data$growth_controls
  nsex <- data$dimensions[["nsex"]]
  nclass <- data$dimensions[["nclass"]]
  max_blocks <- max(controls$n_size_inc_varies)
  out <- array(0, dim = c(nsex, max_blocks, nclass, nclass))

  if (controls$b_use_custom_growth_matrix != 3L) {
    stop("Only BBRKC gamma growth-increment transitions are ported.", call. = FALSE)
  }

  for (sex in seq_len(nsex)) {
    n_size_sex <- data$n_size_sex[sex]
    for (block in seq_len(controls$n_size_inc_varies[sex])) {
      gt <- matrix(0, nrow = nclass, ncol = nclass)
      for (from in seq_len(n_size_sex - 1L)) {
        shape <- growth$molt_increment[sex, block, from] / growth$gscale[sex, block]
        accum <- 0
        for (to in from:(n_size_sex - 1L)) {
          upper_inc <- (data$size_breaks[to + 1L] - data$mid_points[from]) /
            growth$gscale[sex, block]
          cum_inc <- pgamma(as.numeric(upper_inc), shape = as.numeric(shape), scale = 1)
          gt[from, to] <- cum_inc - accum
          accum <- cum_inc
        }
        gt[from, n_size_sex] <- 1 - accum
      }
      gt[n_size_sex, n_size_sex] <- 1
      out[sex, block, , ] <- gt
    }
  }

  out
}

calc_bbrkc_recruitment_size_distribution <- function(model, data) {
  if (!is.null(data$rec_sdd_fixed)) {
    return(data$rec_sdd_fixed)
  }

  controls <- data$growth_controls
  nsex <- data$dimensions[["nsex"]]
  nclass <- data$dimensions[["nclass"]]
  out <- matrix(0, nrow = nsex, ncol = nclass)

  for (sex in seq_len(nsex)) {
    shape <- as.numeric(model$ra[sex] / model$rbeta[sex])
    cdf <- pgamma(data$size_breaks / as.numeric(model$rbeta[sex]), shape = shape, scale = 1)
    probs <- diff(cdf)
    if (controls$max_recruit_size[sex] < nclass) {
      probs[(controls$max_recruit_size[sex] + 1L):nclass] <- 0
    }
    out[sex, ] <- probs / sum(probs)
  }

  out
}

calc_bbrkc_recruits <- function(par, data, model) {
  controls <- data$rec_controls
  years <- data$years
  nsex <- data$dimensions[["nsex"]]
  n_years <- length(years)
  recruits <- matrix(0, nrow = nsex, ncol = n_years)
  rec_dev <- rep(0, n_years)
  logit_rec_prop <- rep(controls$init_logit_sex_ratio, n_years)

  rec_rows <- years >= controls$rdv_syr & years <= controls$rdv_eyr
  rec_dev[rec_rows] <- par$rec_dev_est[seq_len(sum(rec_rows))]
  logit_rec_prop[rec_rows] <- par$logit_rec_prop_est[seq_len(sum(rec_rows))]

  for (iy in seq_len(n_years)) {
    if (controls$initialize_unfished == 1L) {
      total_recruitment <- exp(model$logR0) * nsex
    } else {
      total_recruitment <- exp(model$logRbar) * nsex
    }
    total_recruitment <- total_recruitment * exp(rec_dev[iy])
    if (nsex == 1L) {
      recruits[1, iy] <- total_recruitment
    } else {
      male_prop <- 1 / (1 + exp(-logit_rec_prop[iy]))
      recruits[1, iy] <- total_recruitment * male_prop
      recruits[2, iy] <- total_recruitment * (1 - male_prop)
    }
  }

  list(
    recruits = recruits,
    rec_dev = rec_dev,
    logit_rec_prop = logit_rec_prop,
    rec_sdd = calc_bbrkc_recruitment_size_distribution(model, data)
  )
}

calc_bbrkc_recruitment_residual <- function(par, data, model) {
  n_years <- length(data$years)
  sigR <- exp(par$theta[10])
  sig2R <- 0.5 * sigR * sigR
  out <- vector("list", n_years)
  male_recruits <- model$recruits[1, ]
  out[[1]] <- log(male_recruits[1]) - model$logR0 + sig2R
  base_log_rec <- if (data$rec_controls$initialize_unfished == 1L) {
    model$logR0
  } else {
    model$logRbar
  }
  for (iy in 2:n_years) {
    out[[iy]] <- log(male_recruits[iy]) -
      (1 - model$rho) * base_log_rec -
      model$rho * log(male_recruits[iy - 1L]) +
      sig2R
  }
  do.call(c, out)
}

build_bbrkc_groups <- function(data) {
  dims <- data$dimensions
  grid <- expand.grid(
    shell = seq_len(dims[["nshell"]]),
    maturity = seq_len(dims[["nmature"]]),
    sex = seq_len(dims[["nsex"]]),
    KEEP.OUT.ATTRS = FALSE
  )
  grid <- grid[, c("sex", "maturity", "shell")]
  grid$group <- seq_len(nrow(grid))
  grid
}

calc_bbrkc_growth_block_by_year <- function(data) {
  controls <- data$growth_controls
  years <- data$years
  nsex <- data$dimensions[["nsex"]]
  out <- matrix(1L, nrow = nsex, ncol = length(years))
  for (sex in seq_len(nsex)) {
    if (controls$n_size_inc_varies[sex] <= 1L) {
      next
    }
    for (block in 2:controls$n_size_inc_varies[sex]) {
      out[sex, years >= controls$iyrs_size_changes[[sex]][block - 1L]] <- block
    }
  }
  out
}

calc_bbrkc_initial_numbers <- function(par, data, model) {
  controls <- data$rec_controls
  dims <- data$dimensions
  ngrp <- dims[["nsex"]] * dims[["nmature"]] * dims[["nshell"]]
  nclass <- dims[["nclass"]]

  if (controls$initialize_unfished != 3L) {
    stop("Only FREEPARSSCALED initial numbers are ported for BBRKC.", call. = FALSE)
  }

  logN0 <- matrix(0, nrow = ngrp, ncol = nclass)
  pos <- 13L
  ipnt <- 0L
  for (sex in seq_len(dims[["nsex"]])) {
    for (maturity in seq_len(dims[["nmature"]])) {
      for (shell in seq_len(dims[["nshell"]])) {
        group <- ((sex - 1L) * dims[["nmature"]] * dims[["nshell"]]) +
          ((maturity - 1L) * dims[["nshell"]]) + shell
        for (size_class in seq_len(nclass)) {
          if (ipnt == 0L) {
            logN0[group, size_class] <- 0
          } else {
            logN0[group, size_class] <- par$theta[pos]
            pos <- pos + 1L
          }
          ipnt <- ipnt + 1L
        }
      }
    }
  }

  N0 <- exp(model$logRini + logN0) / sum(exp(logN0))
  list(logN0 = logN0, N0 = N0)
}

calc_bbrkc_population_numbers <- function(par, data, model) {
  dims <- data$dimensions
  nsex <- dims[["nsex"]]
  nmature <- dims[["nmature"]]
  nshell <- dims[["nshell"]]
  nclass <- dims[["nclass"]]
  nseason <- dims[["nseason"]]
  n_years <- length(data$years)
  ngrp <- nsex * nmature * nshell

  if (nshell != 2L || nmature != 1L || data$rec_controls$terminal_molt != 0L) {
    stop("Only BBRKC Case-D population recursion is ported.", call. = FALSE)
  }

  groups <- build_bbrkc_groups(data)
  growth_block <- calc_bbrkc_growth_block_by_year(data)
  initial <- calc_bbrkc_initial_numbers(par, data, model)
  N <- array(0, dim = c(ngrp, n_years + 1L, nseason, nclass))
  N[, 1L, 1L, ] <- initial$N0

  for (iy in seq_len(n_years)) {
    for (season in seq_len(nseason)) {
      next_N <- matrix(0, nrow = ngrp, ncol = nclass)
      last_new_shell_nonmolters <- matrix(0, nrow = nsex, ncol = nclass)
      for (ig in seq_len(ngrp)) {
        sex <- groups$sex[ig]
        maturity <- groups$maturity[ig]
        shell <- groups$shell[ig]
        x <- N[ig, iy, season, ] * model$S[sex, maturity, iy, season, ]

        if (shell == 1L) {
          if (season == data$season_growth) {
            last_new_shell_nonmolters[sex, ] <- x * (1 - model$molt_probability[sex, iy, ])
            x <- as.vector((x * model$molt_probability[sex, iy, ]) %*%
              model$growth_transition[sex, growth_block[sex, iy], , ])
          }
          if (season == data$season_recruitment) {
            x <- x + model$recruits[sex, iy] * model$rec_sdd[sex, ]
          }
          next_N[ig, ] <- next_N[ig, ] + x
        } else {
          z <- rep(0, nclass)
          if (season == data$season_growth) {
            z <- as.vector((x * model$molt_probability[sex, iy, ]) %*%
              model$growth_transition[sex, growth_block[sex, iy], , ])
            x <- x * (1 - model$molt_probability[sex, iy, ]) +
              last_new_shell_nonmolters[sex, ]
          }
          next_N[ig - 1L, ] <- next_N[ig - 1L, ] + z
          next_N[ig, ] <- next_N[ig, ] + x
        }
      }

      if (season == nseason) {
        N[, iy + 1L, 1L, ] <- next_N
      } else {
        N[, iy, season + 1L, ] <- next_N
      }
    }
  }

  list(
    N = N,
    logN0 = initial$logN0,
    N0 = initial$N0,
    groups = groups
  )
}

calc_bbrkc_catch_selectivity <- function(model, fleet, sex, year_index, type, use_discard_report = FALSE) {
  capture <- exp(model$log_slx_capture[fleet, sex, year_index, ])
  retained <- exp(model$log_slx_retaind[fleet, sex, year_index, ])
  switch(
    as.character(type),
    "1" = capture * retained,
    "2" = {
      if (use_discard_report) {
        capture * exp(model$log_slx_discard[fleet, sex, year_index, ])
      } else {
        capture * (1 - retained)
      }
    },
    "0" = capture,
    stop("Unknown catch type: ", type, call. = FALSE)
  )
}

calc_bbrkc_catch_weighted_numbers <- function(model, data, sex, maturity, group, year_index, season, units) {
  nal <- model$N[group, year_index, season, ]
  if (units == 1L) {
    nal <- nal * data$mean_wt[sex, maturity, year_index, ]
  }
  nal
}

calc_bbrkc_catch_mortality_multiplier <- function(model, data, sex, maturity, year_index, season) {
  z <- if (data$season_type[season] == 0L) {
    model$Z2[sex, maturity, year_index, season, ]
  } else {
    model$Z[sex, maturity, year_index, season, ]
  }
  (1 - exp(-z)) / z
}

calc_bbrkc_predicted_catch <- function(data, model) {
  n_series <- length(data$catch$frames)
  zero_ad <- sum(model$N[1, 1, 1, ] * 0)
  catch_rows <- vapply(data$catch$frames, nrow, integer(1))
  groups <- model$groups
  log_q_values <- vector("list", n_series)

  for (series in seq_len(n_series)) {
    frame <- data$catch$frames[[series]]
    q_accumulator <- zero_ad
    nhit <- 0L
    for (row in seq_len(nrow(frame))) {
      cobs <- frame$obs[row]
      effort <- frame$effort[row]
      if (cobs <= 0 || effort <= 0) {
        next
      }
      year_index <- match(frame$year[row], data$years)
      fleet <- as.integer(frame$fleet[row])
      sex_in <- as.integer(frame$sex[row])
      season <- as.integer(frame$season[row])
      type <- as.integer(frame$type[row])
      units <- as.integer(frame$units[row])
      if (sex_in != 0L) {
        q_accumulator <- q_accumulator +
          log(model$ft[fleet, sex_in, year_index, season] / effort)
        nhit <- nhit + 1L
      } else {
        total_obs <- zero_ad
        total_effort <- zero_ad
        for (sex in seq_len(data$dimensions[["nsex"]])) {
          sel <- calc_bbrkc_catch_selectivity(model, fleet, sex, year_index, type)
          for (maturity in seq_len(data$dimensions[["nmature"]])) {
            mult <- calc_bbrkc_catch_mortality_multiplier(model, data, sex, maturity, year_index, season)
            for (group in groups$group[groups$sex == sex & groups$maturity == maturity]) {
              nal <- calc_bbrkc_catch_weighted_numbers(model, data, sex, maturity, group, year_index, season, units)
              total_obs <- total_obs + sum(nal * model$ft[fleet, sex, year_index, season] * sel * mult)
              total_effort <- total_effort + sum(nal * effort * sel * mult)
            }
          }
        }
        q_accumulator <- q_accumulator + log(total_obs / total_effort)
        nhit <- nhit + 1L
      }
    }
    log_q_values[[series]] <- if (nhit > 0L) q_accumulator / nhit else zero_ad
  }
  log_q_catch <- do.call(c, log_q_values)

  pre_values <- vector("list", sum(catch_rows))
  obs_effort_values <- vector("list", sum(catch_rows))
  residual_values <- vector("list", sum(catch_rows))
  out_row <- 1L
  for (series in seq_len(n_series)) {
    frame <- data$catch$frames[[series]]
    for (row in seq_len(nrow(frame))) {
      year_index <- match(frame$year[row], data$years)
      fleet <- as.integer(frame$fleet[row])
      sex_in <- as.integer(frame$sex[row])
      season <- as.integer(frame$season[row])
      type <- as.integer(frame$type[row])
      units <- as.integer(frame$units[row])
      cobs <- frame$obs[row]
      effort <- frame$effort[row]
      sexes <- if (sex_in == 0L) seq_len(data$dimensions[["nsex"]]) else sex_in

      pre <- zero_ad
      obs_effort <- zero_ad
      for (sex in sexes) {
        sel <- calc_bbrkc_catch_selectivity(model, fleet, sex, year_index, type)
        for (maturity in seq_len(data$dimensions[["nmature"]])) {
          mult <- calc_bbrkc_catch_mortality_multiplier(model, data, sex, maturity, year_index, season)
          for (group in groups$group[groups$sex == sex & groups$maturity == maturity]) {
            nal <- calc_bbrkc_catch_weighted_numbers(model, data, sex, maturity, group, year_index, season, units)
            pre <- pre + sum(nal * model$ft[fleet, sex, year_index, season] * sel * mult)
            if (cobs == 0 && effort > 0) {
              obs_effort <- obs_effort +
                sum(nal * exp(log_q_catch[series]) * effort * sel * mult)
            }
          }
        }
      }

      residual <- zero_ad
      if (cobs > 0) {
        residual <- log(cobs) - log(pre)
      } else if (effort > 0) {
        residual <- log(obs_effort) - log(pre)
      }
      pre_values[[out_row]] <- pre
      obs_effort_values[[out_row]] <- obs_effort
      residual_values[[out_row]] <- residual
      out_row <- out_row + 1L
    }
  }

  list(
    pre_catch = do.call(c, pre_values),
    obs_catch_effort = do.call(c, obs_effort_values),
    res_catch = do.call(c, residual_values),
    log_q_catch = log_q_catch
  )
}

calc_bbrkc_predicted_size_compositions <- function(data, model) {
  controls <- data$size_comp_controls
  n_input <- length(data$size_comp$frames)
  pred_input <- vector("list", n_input)
  groups <- model$groups
  zero_ad <- sum(model$N[1, 1, 1, ] * 0)

  for (series in seq_len(n_input)) {
    frame <- data$size_comp$frames[[series]]
    nbin <- data$size_comp$cols[series]
    pred <- vector("list", nrow(frame))
    for (row in seq_len(nrow(frame))) {
      year_index <- match(frame$year[row], data$years)
      if (is.na(year_index) && frame$year[row] == max(data$years) + 1L) {
        year_index <- length(data$years) + 1L
      }
      effect_year_index <- min(year_index, length(data$years))
      fleet <- as.integer(frame$fleet[row])
      season <- as.integer(frame$season[row])
      sex_in <- as.integer(frame$sex[row])
      type <- as.integer(frame$type[row])
      shell_in <- as.integer(frame$shell[row])
      maturity_in <- as.integer(frame$maturity[row])
      sexes <- if (sex_in == 0L) seq_len(data$dimensions[["nsex"]]) else sex_in
      dntmp <- model$N[1, year_index, season, ] * 0

      for (sex in sexes) {
        sel <- calc_bbrkc_catch_selectivity(
          model, fleet, sex, effect_year_index, type,
          use_discard_report = TRUE
        )
        for (maturity in seq_len(data$dimensions[["nmature"]])) {
          if (maturity_in != 0L && maturity_in != maturity) {
            next
          }
          mult <- if (controls$lf_catch_in[series] == 1L) {
            calc_bbrkc_catch_mortality_multiplier(model, data, sex, maturity, effect_year_index, season)
          } else {
            rep(1, data$dimensions[["nclass"]])
          }
          for (group in groups$group[groups$sex == sex & groups$maturity == maturity]) {
            shell <- groups$shell[group]
            if (shell_in != 0L && shell_in != shell) {
              next
            }
            dntmp <- dntmp + model$N[group, year_index, season, ] * sel * mult
          }
        }
      }
      pred[[row]] <- if (nbin < data$dimensions[["nclass"]]) {
        c(dntmp[seq_len(nbin - 1L)], sum(dntmp[nbin:data$dimensions[["nclass"]]]))
      } else {
        dntmp
      }
    }
    pred_input[[series]] <- pred
  }

  n_modified <- max(controls$comp_aggregator)
  pred_modified <- vector("list", n_modified)
  comp_rows <- integer(n_modified)
  comp_cols <- integer(n_modified)
  for (modified in seq_len(n_modified)) {
    input_series <- which(controls$comp_aggregator == modified)
    nrow_modified <- nrow(data$size_comp$frames[[input_series[1]]])
    comp_rows[modified] <- nrow_modified
    comp_cols[modified] <- sum(data$size_comp$cols[input_series])
    pred_rows <- vector("list", nrow_modified)
    for (row in seq_len(nrow_modified)) {
      row_values <- vector("list", length(input_series))
      for (i in seq_along(input_series)) {
        row_values[[i]] <- pred_input[[input_series[i]]][[row]]
      }
      row_pred <- do.call(c, row_values)
      row_sum <- sum(row_pred) + zero_ad
      pred_rows[[row]] <- row_pred / row_sum
    }
    pred_modified[[modified]] <- do.call(c, pred_rows)
  }

  list(
    pre_size_comps = do.call(c, pred_modified),
    size_comp_rows = comp_rows,
    size_comp_cols = comp_cols
  )
}

calc_bbrkc_survey_weighted_numbers <- function(model, data, sex, maturity, group, year_index, season, units) {
  effect_year_index <- min(year_index, length(data$years))
  numbers <- model$N[group, year_index, season, ]
  if (units == 1L) {
    numbers * data$mean_wt[sex, maturity, effect_year_index, ]
  } else {
    numbers
  }
}

calc_bbrkc_predicted_survey <- function(par, data, model) {
  survey <- data$survey$data
  controls <- data$q_controls
  groups <- model$groups
  n_rows <- nrow(survey)
  n_surveys <- length(data$survey$index_type)
  zero_ad <- sum(model$N[1, 1, 1, ] * 0)
  vulnerable <- vector("list", n_rows)
  pre_cpue <- vector("list", n_rows)
  res_cpue <- vector("list", n_rows)
  survey_q <- vector("list", n_surveys)

  for (row in seq_len(n_rows)) {
    year_index <- match(survey$year[row], data$years)
    if (is.na(year_index) && survey$year[row] == max(data$years) + 1L) {
      year_index <- length(data$years) + 1L
    }
    effect_year_index <- min(year_index, length(data$years))
    season <- as.integer(survey$season[row])
    fleet <- as.integer(survey$fleet[row])
    sex_in <- as.integer(survey$sex[row])
    maturity_in <- as.integer(survey$maturity[row])
    units <- as.integer(survey$units[row])
    cpue_time <- survey$cpue_time[row]
    series <- as.integer(survey$index[row])
    sexes <- if (sex_in == 0L) seq_len(data$dimensions[["nsex"]]) else sex_in
    row_v <- zero_ad

    for (sex in sexes) {
      sel <- exp(model$log_slx_capture[fleet, sex, effect_year_index, ])
      ret <- exp(model$log_slx_retaind[fleet, sex, effect_year_index, ])
      survey_availability <- if (data$survey$index_type[series] == 2L) sel * ret else sel
      if (cpue_time > 0) {
        z_at_time <- if (data$season_type[season] == 1L) {
          model$Z2[sex, 1, effect_year_index, season, ] * cpue_time
        } else {
          model$Z[sex, 1, effect_year_index, season, ] * cpue_time + 1e-10
        }
      } else {
        z_at_time <- rep(0, data$dimensions[["nclass"]])
      }

      for (maturity in seq_len(data$dimensions[["nmature"]])) {
        if (maturity_in != 0L && maturity_in != maturity) {
          next
        }
        for (group in groups$group[groups$sex == sex & groups$maturity == maturity]) {
          nal <- calc_bbrkc_survey_weighted_numbers(
            model, data, sex, maturity, group, year_index, season, units
          )
          row_v <- row_v + sum(nal * survey_availability * exp(-z_at_time))
        }
      }
    }
    vulnerable[[row]] <- row_v
  }

  vulnerable_vec <- do.call(c, vulnerable)
  for (series in seq_len(n_surveys)) {
    rows <- which(survey$index == series)
    if (controls$q_anal[series] == 1L) {
      ztot1 <- zero_ad
      ztot2 <- 0
      if (controls$prior_qtype[series] == 2L) {
        ztot1 <- ztot1 + log(controls$prior_p1[series]) / controls$prior_p2[series]^2
        ztot2 <- ztot2 + 1 / controls$prior_p2[series]^2
      }
      for (row in rows) {
        cvadd2 <- if (controls$add_cv_links[series] > 0L) {
          log(1 + exp(par$log_add_cv[controls$add_cv_links[series]])^2)
        } else {
          0
        }
        cvobs2 <- log(1 + survey$cv[row]^2) / controls$cpue_lambda[series]
        variance <- cvobs2 + cvadd2
        zt <- log(survey$obs[row]) - log(vulnerable_vec[row])
        ztot1 <- ztot1 + zt / variance
        ztot2 <- ztot2 + 1 / variance
      }
      survey_q[[series]] <- exp(ztot1 / ztot2)
    } else {
      survey_q[[series]] <- par$survey_q[series]
    }
  }
  survey_q_vec <- do.call(c, survey_q)

  for (row in seq_len(n_rows)) {
    series <- as.integer(survey$index[row])
    pred <- survey_q_vec[series] * vulnerable_vec[row]
    pre_cpue[[row]] <- pred
    res_cpue[[row]] <- log(survey$obs[row]) - log(pred)
  }

  list(
    vulnerable_cpue = vulnerable_vec,
    pre_cpue = do.call(c, pre_cpue),
    res_cpue = do.call(c, res_cpue),
    survey_q_calc = survey_q_vec
  )
}

initialize_bbrkc_model_parameters <- function(par, data, predictions = TRUE) {
  theta <- par$theta
  m0 <- c(theta[1], theta[1] * exp(theta[2]))
  retention <- calc_bbrkc_retention_selectivity(par$log_slx_pars, par$Asymret, data)
  model <- list(
    M0 = m0,
    logR0 = theta[3],
    logRini = theta[4],
    logRbar = theta[5],
    ra = c(theta[6], theta[6] * exp(theta[8])),
    rbeta = c(theta[7], theta[7] * exp(theta[9])),
    logSigmaR = theta[10],
    steepness = theta[11],
    rho = theta[12],
    log_slx_capture = calc_bbrkc_capture_selectivity(par$log_slx_pars, data),
    log_slx_retaind = retention$log_slx_retaind,
    log_slx_discard = retention$log_slx_discard
  )
  model$M <- calc_bbrkc_natural_mortality(par, data, model)
  fishing <- calc_bbrkc_fishing_mortality(par, data, model)
  model <- c(model, fishing)
  total_mortality <- calc_bbrkc_total_mortality(data, model)
  model <- c(model, total_mortality)
  growth_par <- par
  if (!is.null(data$growth_parameters_initial)) {
    growth_par$Grwth <- data$growth_parameters_initial
  }
  growth <- unpack_bbrkc_growth_parameters(growth_par, data)
  growth$molt_probability <- calc_bbrkc_molting_probability(growth, data)
  growth$growth_transition <- calc_bbrkc_growth_transition(growth, data)
  model <- c(model, growth)
  recruitment <- calc_bbrkc_recruits(par, data, model)
  model <- c(model, recruitment)
  model$res_recruit <- calc_bbrkc_recruitment_residual(par, data, model)
  population <- calc_bbrkc_population_numbers(par, data, model)
  model <- c(model, population)
  if (predictions) {
    catch_prediction <- calc_bbrkc_predicted_catch(data, model)
    model <- c(model, catch_prediction)
    size_prediction <- calc_bbrkc_predicted_size_compositions(data, model)
    model <- c(model, size_prediction)
    survey_prediction <- calc_bbrkc_predicted_survey(par, data, model)
    model <- c(model, survey_prediction)
  }
  model
}

multinomial_nll <- function(obs, pred, log_effn) {
  eps <- 1e-12
  nll <- 0
  for (i in seq_len(nrow(obs))) {
    p <- pred[i, ] / sum(pred[i, ])
    o <- obs[i, ] / sum(obs[i, ])
    vn <- exp(log_effn[i])
    sobs <- vn * o
    nll <- nll - lgamma(vn)
    for (j in seq_along(o)) {
      if (o[j] > 0) {
        nll <- nll + lgamma(sobs[j])
      }
    }
    nll <- nll - sum(sobs * log(eps + p))
  }
  nll
}

multinomial_nll_flat <- function(obs, pred, log_effn) {
  eps <- 1e-12
  nll <- 0
  cols <- ncol(obs)
  for (i in seq_len(nrow(obs))) {
    row_start <- (i - 1L) * cols
    psum <- 0
    for (j in seq_len(cols)) {
      psum <- psum + pred[row_start + j]
    }
    o <- obs[i, ] / sum(obs[i, ])
    vn <- exp(log_effn[i])
    nll <- nll - lgamma(vn)
    for (j in seq_len(cols)) {
      sobs <- vn * o[j]
      if (o[j] > 0) {
        nll <- nll + lgamma(sobs)
      }
      nll <- nll - sobs * log(eps + pred[row_start + j] / psum)
    }
  }
  nll
}

robust_multinomial_nll <- function(obs, pred, log_effn) {
  eps <- 1e-08
  nll <- 0
  a <- 0.1 / ncol(obs)
  effn <- exp(log_effn)
  for (i in seq_len(nrow(obs))) {
    o <- obs[i, ] + eps
    o <- o / sum(o)
    psum <- sum(pred[i, ] + eps)
    for (j in seq_along(o)) {
      p <- (pred[i, j] + eps) / psum
      v <- a + o[j] * (1 - o[j])
      l <- 0.5 * ((p - o[j])^2 / v)
      nll <- nll - log(exp(-effn[i] * l) + 0.01)
      nll <- nll + 0.5 * log(v / effn[i])
    }
  }
  nll
}

robust_multinomial_nll_flat <- function(obs, pred, log_effn) {
  eps <- 1e-08
  nll <- 0
  cols <- ncol(obs)
  a <- 0.1 / cols
  effn <- exp(log_effn)
  for (i in seq_len(nrow(obs))) {
    row_start <- (i - 1L) * cols
    psum <- 0
    for (j in seq_len(cols)) {
      psum <- psum + pred[row_start + j] + eps
    }
    o <- obs[i, ] + eps
    o <- o / sum(o)
    for (j in seq_len(cols)) {
      p <- (pred[row_start + j] + eps) / psum
      v <- a + o[j] * (1 - o[j])
      l <- 0.5 * ((p - o[j])^2 / v)
      nll <- nll - log(exp(-effn[i] * l) + 0.01)
      nll <- nll + 0.5 * log(v / effn[i])
    }
  }
  nll
}

dirichlet_nll <- function(obs, pred, log_effn) {
  eps <- 1e-10
  nll <- 0
  for (i in seq_len(nrow(obs))) {
    p <- pred[i, ] / sum(pred[i, ])
    o <- obs[i, ] / sum(obs[i, ])
    alpha <- exp(log_effn[i]) * p
    nll <- nll - (sum((alpha - 1) * log(eps + o)) -
      (sum(lgamma(alpha)) - lgamma(sum(alpha))))
  }
  nll
}

dirichlet_nll_flat <- function(obs, pred, log_effn) {
  eps <- 1e-10
  nll <- 0
  cols <- ncol(obs)
  for (i in seq_len(nrow(obs))) {
    row_start <- (i - 1L) * cols
    psum <- 0
    for (j in seq_len(cols)) {
      psum <- psum + pred[row_start + j]
    }
    o <- obs[i, ] / sum(obs[i, ])
    alpha0 <- 0
    lmnB <- 0
    sj <- 0
    for (j in seq_len(cols)) {
      alpha <- exp(log_effn[i]) * pred[row_start + j] / psum
      alpha0 <- alpha0 + alpha
      lmnB <- lmnB + lgamma(alpha)
      sj <- sj + (alpha - 1) * log(eps + o[j])
    }
    lmnB <- lmnB - lgamma(alpha0)
    nll <- nll - (sj - lmnB)
  }
  nll
}

build_bbrkc_observed_size_compositions <- function(data) {
  controls <- data$size_comp_controls
  n_modified <- max(controls$comp_aggregator)
  observed <- vector("list", n_modified)
  sample_size <- vector("list", n_modified)
  likelihood_type <- integer(n_modified)
  lambda <- numeric(n_modified)
  emphasis <- numeric(n_modified)
  tail_compression <- logical(n_modified)

  for (modified in seq_len(n_modified)) {
    input_series <- which(controls$comp_aggregator == modified)
    nrow_modified <- nrow(data$size_comp$frames[[input_series[1]]])
    observed_rows <- vector("list", nrow_modified)
    sample_size[[modified]] <- rep(0, nrow_modified)
    for (row in seq_len(nrow_modified)) {
      row_values <- vector("list", length(input_series))
      for (i in seq_along(input_series)) {
        series <- input_series[i]
        frame <- data$size_comp$frames[[series]]
        bins <- grep("^bin_", names(frame))
        row_values[[i]] <- as.numeric(frame[row, bins])
        sample_size[[modified]][row] <- sample_size[[modified]][row] + frame$nsamp[row]
      }
      observed_rows[[row]] <- do.call(c, row_values)
    }
    observed[[modified]] <- do.call(rbind, observed_rows)
    likelihood_type[modified] <- controls$likelihood_type_in[input_series[1]]
    lambda[modified] <- controls$lambda[input_series[1]]
    emphasis[modified] <- controls$emphasis[input_series[1]]
    tail_compression[modified] <- controls$tail_compression_in[input_series[1]] > 0
  }

  list(
    observed = observed,
    sample_size = sample_size,
    likelihood_type = likelihood_type,
    lambda = lambda,
    emphasis = emphasis,
    tail_compression = tail_compression
  )
}

calc_bbrkc_penalties <- function(par, data, model, zero_ad) {
  out <- vector("list", 13)
  for (i in seq_along(out)) {
    out[[i]] <- zero_ad
  }

  log_fdev <- split_parameter_by_fleet(par$log_fdev, data$fdev_indices)
  log_fdov <- split_parameter_by_fleet(par$log_fdov, data$fdov_indices)
  for (fleet in seq_along(log_fdev)) {
    if (length(log_fdev[[fleet]]) > 0L) {
      out[[1]] <- out[[1]] +
        data$likelihood_emphasis$fdev_penalty[fleet, 1] * mean(log_fdev[[fleet]])^2
      out[[11]] <- out[[11]] +
        data$likelihood_emphasis$fdev_penalty[fleet, 3] * sum(log_fdev[[fleet]]^2)
    }
    if (length(log_fdov[[fleet]]) > 0L) {
      out[[1]] <- out[[1]] +
        data$likelihood_emphasis$fdev_penalty[fleet, 2] * mean(log_fdov[[fleet]])^2
      out[[12]] <- out[[12]] +
        data$likelihood_emphasis$fdev_penalty[fleet, 4] * sum(log_fdov[[fleet]]^2)
    }
  }

  if (length(par$m_dev_est) > 0L) {
    out[[3]] <- admb_dnorm_nll(par$m_dev_est, data$m_controls$mdev_sd)
  }
  if (length(par$rec_dev_est) > 0L) {
    out[[6]] <- first_difference_nll(model$rec_dev, 1)
  }
  if (data$dimensions[["nsex"]] > 1L) {
    out[[7]] <- (log(sum(model$recruits[2, ])) - log(sum(model$recruits[1, ])))^2
  }
  for (group in seq_len(nrow(model$logN0))) {
    out[[10]] <- out[[10]] + first_difference_nll(model$logN0[group, ], 1)
  }

  out <- do.call(c, out)
  names(out) <- paste0("penalty_", seq_along(out))
  out
}

calc_bbrkc_prior_density <- function(par, data, zero_ad) {
  out <- rep(0, length(data$admb_reference$prior_density))
  iprior <- 1L

  set_prior <- function(value) {
    out[iprior] <<- value
    iprior <<- iprior + 1L
  }
  skip_prior <- function(n) {
    if (n > 0L) {
      iprior <<- iprior + n
    }
  }
  active_count <- function(x) sum(x > 0)

  controls <- data$prior_controls

  for (i in seq_along(par$theta)) {
    row <- controls$theta[i, ]
    if (row[4] > 0) {
      x <- par$theta[i]
      if (row[5] == 3) {
        x <- (x - row[2]) / (row[3] - row[2])
      }
      set_prior(admb_prior_nll(row[5], x, row[6], row[7]))
    }
  }

  for (i in seq_along(par$Grwth)) {
    row <- controls$Grwth[i, ]
    if (row[4] > 0) {
      x <- par$Grwth[i]
      if (row[5] == 3) {
        x <- (x - row[2]) / (row[3] - row[2])
      }
      set_prior(admb_prior_nll(row[5], x, row[6], row[7]))
    }
  }

  for (i in seq_along(par$log_slx_pars)) {
    row <- controls$log_slx_pars[i, ]
    if (row[11] > 0) {
      x <- exp(par$log_slx_pars[i])
      p1 <- if (row[8] == 0) row[6] else row[9]
      p2 <- if (row[8] == 0) row[7] else row[10]
      set_prior(admb_prior_nll(row[8], x, p1, p2))
    }
  }

  for (i in seq_along(par$Asymret)) {
    row <- controls$Asymret[i, ]
    if (row[7] > 0) {
      set_prior(admb_prior_nll(0, par$Asymret[i], row[5], row[6]))
    }
  }

  skip_prior(active_count(data$parameter_phases$log_fbar))
  skip_prior(length(unlist(data$fdev_indices[data$parameter_phases$log_fbar > 0], use.names = FALSE)))
  skip_prior(active_count(data$parameter_phases$log_foff))
  skip_prior(length(unlist(data$fdov_indices[data$parameter_phases$log_foff > 0], use.names = FALSE)))

  if (data$rec_controls$rec_ini_phz > 0 && data$rec_controls$rdv_phz > 0) {
    skip_prior(length(par$rec_ini))
  }
  if (data$rec_controls$rdv_phz > 0) {
    skip_prior(length(par$rec_dev_est))
  }
  if (data$rec_controls$rec_prop_phz > 0) {
    skip_prior(length(par$logit_rec_prop_est))
  }

  for (i in seq_along(par$m_dev_est)) {
    row <- controls$m_dev_est[i, ]
    if (row[4] > 0 && row[5] >= 0) {
      set_prior(admb_prior_nll(0, par$m_dev_est[i], row[2], row[3]))
    }
  }

  skip_prior(active_count(data$size_comp_controls$nvn_phz_in))

  for (i in seq_along(par$survey_q)) {
    row <- controls$q[i, ]
    if (row[4] > 0) {
      value <- admb_prior_nll(row[5], par$survey_q[i], row[6], row[7])
      set_prior(value)
    }
  }

  for (i in seq_along(par$log_add_cv)) {
    row <- controls$add_cv[i, ]
    if (row[4] > 0) {
      p1 <- if (row[5] == 0) row[2] else row[6]
      p2 <- if (row[5] == 0) row[3] else row[7]
      set_prior(admb_prior_nll(row[5], exp(par$log_add_cv[i]), p1, p2))
    }
  }

  for (i in seq_along(par$m_mat_mult)) {
    row <- controls$m_mat_mult[i, ]
    if (row[4] > 0) {
      p1 <- if (row[5] == 0) row[2] else row[6]
      p2 <- if (row[5] == 0) row[3] else row[7]
      set_prior(admb_prior_nll(row[5], par$m_mat_mult[i], p1, p2))
    }
  }

  out + zero_ad
}

make_gmacs_rtmb_data <- function(inputs) {
  dat <- inputs$data
  dims <- dat$dimensions
  fdev_labels <- attr(inputs$parameters$log_fdev, "parameter_names")
  fdov_labels <- attr(inputs$parameters$log_fdov, "parameter_names")
  data <- list(
    dimensions = dims,
    years = seq.int(dims[["syr"]], dims[["nyr"]]),
    fleet_names = dat$fleet_names,
    fdev_indices = lapply(dat$fleet_names, function(fleet) {
      which(startsWith(fdev_labels, paste0("Log_fdev_", fleet)))
    }),
    fdov_indices = lapply(dat$fleet_names, function(fleet) {
      which(startsWith(fdov_labels, paste0("Log_fdov_", fleet)))
    }),
    mid_points = dat$mid_points,
    size_breaks = dat$size_breaks,
    n_size_sex = dat$n_size_sex,
    mean_wt = {
      years <- seq.int(dims[["syr"]], dims[["nyr"]])
      x <- array(0, dim = c(dims[["nsex"]], dims[["nmature"]], length(years), dims[["nclass"]]))
      for (sex in seq_len(dims[["nsex"]])) {
        for (maturity in seq_len(dims[["nmature"]])) {
          for (iy in seq_along(years)) {
            x[sex, maturity, iy, ] <- dat$allometry$mean_wt[sex, maturity, ]
          }
        }
      }
      x
    },
    maturity = dat$allometry$maturity,
    legal = dat$allometry$legal,
    season_recruitment = dat$season_recruitment,
    season_growth = dat$season_growth,
    season_ssb = dat$season_ssb,
    season_N = dat$season_N,
    m_prop = dat$m_prop,
    m_controls = inputs$control_phases$m_controls,
    growth_controls = inputs$control_phases$growth_controls,
    size_comp_controls = inputs$control_phases$size_comp_controls,
    rec_controls = inputs$control_phases$rec_controls,
    q_controls = inputs$control_phases$q_controls,
    parameter_phases = inputs$control_phases,
    prior_controls = inputs$control_phases$prior_controls,
    likelihood_emphasis = inputs$control_phases$likelihood_emphasis,
    season_type = dat$season_type,
    catch = dat$catch,
    survey = dat$survey,
    size_comp = dat$size_comp,
    growth = dat$growth,
    environment = dat$environment,
    growth_parameters_initial = {
      x <- inputs$parameters$Grwth
      attributes(x) <- NULL
      x
    },
    admb_reference = inputs$admb_reference
  )
  initial_parameters <- make_gmacs_parameter_list(inputs)
  initial_growth <- unpack_bbrkc_growth_parameters(initial_parameters, data)
  data$growth_transition_fixed <- calc_bbrkc_growth_transition(initial_growth, data)
  theta <- initial_parameters$theta
  initial_model <- list(
    ra = c(theta[6], theta[6] * exp(theta[8])),
    rbeta = c(theta[7], theta[7] * exp(theta[9]))
  )
  data$rec_sdd_fixed <- calc_bbrkc_recruitment_size_distribution(initial_model, data)
  data$observed_size_compositions <- build_bbrkc_observed_size_compositions(data)
  data
}

make_gmacs_parameter_list <- function(inputs) {
  lapply(inputs$parameters, function(x) {
    attributes(x) <- NULL
    x
  })
}

active_factor <- function(active) {
  out <- rep(NA_integer_, length(active))
  out[active] <- seq_len(sum(active))
  factor(out)
}

make_bbrkc_rtmb_map <- function(inputs) {
  pars <- inputs$parameters
  phases <- inputs$control_phases
  fleets <- inputs$data$fleet_names

  map <- list(
    theta = active_factor(phases$theta > 0),
    Grwth = active_factor(phases$Grwth > 0),
    log_slx_pars = active_factor(phases$log_slx_pars > 0),
    Asymret = active_factor(phases$Asymret > 0),
    slx_env_pars = active_factor(rep(FALSE, length(pars$slx_env_pars))),
    sel_devs = active_factor(rep(FALSE, length(pars$sel_devs))),
    log_fbar = active_factor(phases$log_fbar > 0),
    log_foff = active_factor(phases$log_foff > 0),
    m_dev_est = active_factor(phases$m_dev_est > 0),
    m_mat_mult = active_factor(phases$m_mat_mult > 0),
    log_vn = active_factor(phases$log_vn > 0),
    survey_q = active_factor(phases$survey_q > 0),
    log_add_cv = active_factor(phases$log_add_cv > 0),
    rec_ini = active_factor(rep(phases$rec_ini_phz > 0, length(pars$rec_ini))),
    rec_dev_est = active_factor(rep(phases$rdv_phz > 0, length(pars$rec_dev_est)))
  )

  fdev_labels <- attr(pars$log_fdev, "parameter_names")
  fdov_labels <- attr(pars$log_fdov, "parameter_names")
  fdev_active <- rep(FALSE, length(pars$log_fdev))
  fdov_active <- rep(FALSE, length(pars$log_fdov))
  for (i in seq_along(fleets)) {
    fdev_active <- fdev_active |
      (startsWith(fdev_labels, paste0("Log_fdev_", fleets[i])) & phases$log_fbar[i] > 0)
    fdov_active <- fdov_active |
      (startsWith(fdov_labels, paste0("Log_fdov_", fleets[i])) & phases$log_foff[i] > 0)
  }
  map$log_fdev <- active_factor(fdev_active)
  map$log_fdov <- active_factor(fdov_active)

  logit_active <- rep(phases$rec_prop_phz > 0, length(pars$logit_rec_prop_est))
  if (any(logit_active)) {
    logit_active[length(logit_active)] <- FALSE
  }
  map$logit_rec_prop_est <- active_factor(logit_active)

  map[names(pars)]
}

count_active_map_parameters <- function(map) {
  sum(vapply(map, function(x) length(levels(x)), integer(1)))
}

selectivity_array_to_data_frame <- function(x, data) {
  out_array <- exp(x)
  grid <- expand.grid(
    fleet = data$fleet_names,
    sex = c("male", "female"),
    year = data$years,
    size_class = seq_along(data$mid_points),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$value <- as.vector(out_array)
  grid
}

compare_selectivity_to_admb <- function(model, inputs, component) {
  data <- make_gmacs_rtmb_data(inputs)
  rtmb_array <- switch(
    component,
    capture = model$log_slx_capture,
    retained = model$log_slx_retaind,
    discard = model$log_slx_discard,
    stop("Unknown selectivity component: ", component, call. = FALSE)
  )
  rtmb <- selectivity_array_to_data_frame(rtmb_array, data)
  admb <- inputs$admb_selectivity[[component]]
  merged <- merge(
    admb,
    rtmb,
    by = c("year", "sex", "fleet", "size_class"),
    suffixes = c("_admb", "_rtmb")
  )
  merged$abs_diff <- abs(merged$value_admb - merged$value_rtmb)
  list(
    max_abs_diff = max(merged$abs_diff),
    mean_abs_diff = mean(merged$abs_diff),
    n = nrow(merged),
    data = merged
  )
}

fully_selected_f_to_data_frame <- function(ft, data) {
  grid <- expand.grid(
    fleet = seq_len(data$dimensions[["nfleet"]]),
    sex = c("male", "female"),
    year = data$years,
    season = seq_len(data$dimensions[["nseason"]]),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$value <- as.vector(ft)
  grid
}

compare_fully_selected_f_to_admb <- function(model, inputs) {
  data <- make_gmacs_rtmb_data(inputs)
  rtmb <- fully_selected_f_to_data_frame(model$ft, data)
  admb <- inputs$admb_fishing_mortality$fully_selected_fleet
  merged <- merge(
    admb,
    rtmb,
    by = c("sex", "year", "season", "fleet"),
    suffixes = c("_admb", "_rtmb")
  )
  merged$abs_diff <- abs(merged$value_admb - merged$value_rtmb)
  list(
    max_abs_diff = max(merged$abs_diff),
    mean_abs_diff = mean(merged$abs_diff),
    n = nrow(merged),
    data = merged
  )
}

f_at_size_to_data_frame <- function(F, data) {
  grid <- expand.grid(
    sex = c("male", "female"),
    year = data$years,
    season = seq_len(data$dimensions[["nseason"]]),
    size_class = seq_len(data$dimensions[["nclass"]]),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$value <- as.vector(F)
  grid
}

compare_f_at_size_to_admb <- function(model, inputs) {
  data <- make_gmacs_rtmb_data(inputs)
  rtmb <- f_at_size_to_data_frame(model$F, data)
  admb <- inputs$admb_fishing_mortality$at_size_continuous
  merged <- merge(
    admb,
    rtmb,
    by = c("sex", "year", "season", "size_class"),
    suffixes = c("_admb", "_rtmb")
  )
  merged$abs_diff <- abs(merged$value_admb - merged$value_rtmb)
  list(
    max_abs_diff = max(merged$abs_diff),
    mean_abs_diff = mean(merged$abs_diff),
    n = nrow(merged),
    data = merged
  )
}

natural_mortality_to_data_frame <- function(M, data) {
  grid <- expand.grid(
    sex = c("male", "female"),
    maturity = seq_len(data$dimensions[["nmature"]]),
    year = data$years,
    size_class = seq_len(data$dimensions[["nclass"]]),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$value <- as.vector(M)
  grid
}

compare_natural_mortality_to_admb <- function(model, inputs) {
  data <- make_gmacs_rtmb_data(inputs)
  rtmb <- natural_mortality_to_data_frame(model$M, data)
  admb <- inputs$admb_natural_mortality
  merged <- merge(
    admb,
    rtmb,
    by = c("sex", "maturity", "year", "size_class"),
    suffixes = c("_admb", "_rtmb")
  )
  merged$abs_diff <- abs(merged$value_admb - merged$value_rtmb)
  list(
    max_abs_diff = max(merged$abs_diff),
    mean_abs_diff = mean(merged$abs_diff),
    n = nrow(merged),
    data = merged
  )
}

total_mortality_to_data_frame <- function(Z, data) {
  grid <- expand.grid(
    sex = c("male", "female"),
    maturity = seq_len(data$dimensions[["nmature"]]),
    year = data$years,
    season = seq_len(data$dimensions[["nseason"]]),
    size_class = seq_len(data$dimensions[["nclass"]]),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$value <- as.vector(Z)
  grid
}

compare_total_mortality_to_admb <- function(model, inputs, component) {
  data <- make_gmacs_rtmb_data(inputs)
  rtmb_array <- switch(
    component,
    continuous = model$Z,
    discrete = model$Z2,
    stop("Unknown total mortality component: ", component, call. = FALSE)
  )
  rtmb <- total_mortality_to_data_frame(rtmb_array, data)
  admb <- inputs$admb_total_mortality[[component]]
  merged <- merge(
    admb,
    rtmb,
    by = c("sex", "maturity", "year", "season", "size_class"),
    suffixes = c("_admb", "_rtmb")
  )
  merged$abs_diff <- abs(merged$value_admb - merged$value_rtmb)
  list(
    max_abs_diff = max(merged$abs_diff),
    mean_abs_diff = mean(merged$abs_diff),
    n = nrow(merged),
    data = merged
  )
}

molt_probability_to_data_frame <- function(x, data) {
  grid <- expand.grid(
    sex = c("male", "female"),
    year = data$years,
    size_class = seq_len(data$dimensions[["nclass"]]),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$value <- as.vector(x)
  grid
}

compare_molt_probability_to_admb <- function(model, inputs) {
  data <- make_gmacs_rtmb_data(inputs)
  rtmb <- molt_probability_to_data_frame(model$molt_probability, data)
  admb <- inputs$admb_growth$molt_probability
  merged <- merge(
    admb,
    rtmb,
    by = c("sex", "year", "size_class"),
    suffixes = c("_admb", "_rtmb")
  )
  merged$abs_diff <- abs(merged$value_admb - merged$value_rtmb)
  list(
    max_abs_diff = max(merged$abs_diff),
    mean_abs_diff = mean(merged$abs_diff),
    n = nrow(merged),
    data = merged
  )
}

growth_transition_to_data_frame <- function(x, data) {
  nsex <- data$dimensions[["nsex"]]
  nclass <- data$dimensions[["nclass"]]
  records <- vector("list", 0)
  k <- 1L
  for (sex in seq_len(nsex)) {
    sex_name <- c("male", "female")[sex]
    for (increment_no in seq_len(data$growth_controls$n_size_inc_varies[sex])) {
      mat <- x[sex, increment_no, , ]
      records[[k]] <- data.frame(
        sex = sex_name,
        increment_no = increment_no,
        to_size_class = rep(seq_len(nclass), each = nclass),
        from_size_class = rep(seq_len(nclass), times = nclass),
        value = as.vector(mat),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, records)
}

compare_growth_transition_to_admb <- function(model, inputs) {
  data <- make_gmacs_rtmb_data(inputs)
  rtmb <- growth_transition_to_data_frame(model$growth_transition, data)
  admb <- inputs$admb_growth$growth_transition
  merged <- merge(
    admb,
    rtmb,
    by = c("sex", "increment_no", "to_size_class", "from_size_class"),
    suffixes = c("_admb", "_rtmb")
  )
  merged$abs_diff <- abs(merged$value_admb - merged$value_rtmb)
  list(
    max_abs_diff = max(merged$abs_diff),
    mean_abs_diff = mean(merged$abs_diff),
    n = nrow(merged),
    data = merged
  )
}

recruits_to_data_frame <- function(x, data) {
  grid <- expand.grid(
    sex = c("male", "female"),
    year = data$years,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$value <- as.vector(x)
  grid
}

compare_recruits_to_admb <- function(model, inputs) {
  data <- make_gmacs_rtmb_data(inputs)
  rtmb <- recruits_to_data_frame(model$recruits, data)
  admb_summary <- inputs$admb_overall_summary
  admb <- rbind(
    data.frame(
      sex = "male",
      year = admb_summary$Year,
      value = admb_summary$Recruit_male,
      stringsAsFactors = FALSE
    ),
    data.frame(
      sex = "female",
      year = admb_summary$Year,
      value = admb_summary$Recruit_female,
      stringsAsFactors = FALSE
    )
  )
  merged <- merge(
    admb,
    rtmb,
    by = c("sex", "year"),
    suffixes = c("_admb", "_rtmb")
  )
  merged$abs_diff <- abs(merged$value_admb - merged$value_rtmb)
  list(
    max_abs_diff = max(merged$abs_diff),
    mean_abs_diff = mean(merged$abs_diff),
    n = nrow(merged),
    data = merged
  )
}

rec_sdd_to_data_frame <- function(x, data) {
  grid <- expand.grid(
    sex = c("male", "female"),
    size_class = seq_len(data$dimensions[["nclass"]]),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$value <- as.vector(x)
  grid
}

numbers_at_size_to_data_frame <- function(model, data, component) {
  season <- data$season_N
  groups <- model$groups
  N <- model$N[, seq_along(data$years), season, , drop = FALSE]
  nclass <- data$dimensions[["nclass"]]

  group_index <- switch(
    component,
    total = seq_len(nrow(groups)),
    males = groups$group[groups$sex == 1L],
    females = groups$group[groups$sex == 2L],
    males_new = groups$group[groups$sex == 1L & groups$shell == 1L],
    females_new = groups$group[groups$sex == 2L & groups$shell == 1L],
    males_old = groups$group[groups$sex == 1L & groups$shell == 2L],
    females_old = groups$group[groups$sex == 2L & groups$shell == 2L],
    stop("Unknown numbers-at-size component: ", component, call. = FALSE)
  )

  records <- vector("list", length(data$years))
  for (iy in seq_along(data$years)) {
    values <- as.vector(apply(N[group_index, iy, 1L, , drop = FALSE], 4L, sum))
    records[[iy]] <- data.frame(
      year = data$years[iy],
      size_class = seq_len(nclass),
      value = values,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, records)
}

compare_numbers_at_size_to_admb <- function(model, inputs, component) {
  data <- make_gmacs_rtmb_data(inputs)
  rtmb <- numbers_at_size_to_data_frame(model, data, component)
  admb <- inputs$admb_numbers_at_size[[component]]
  merged <- merge(
    admb,
    rtmb,
    by = c("year", "size_class"),
    suffixes = c("_admb", "_rtmb")
  )
  merged$abs_diff <- abs(merged$value_admb - merged$value_rtmb)
  list(
    max_abs_diff = max(merged$abs_diff),
    mean_abs_diff = mean(merged$abs_diff),
    n = nrow(merged),
    data = merged
  )
}

catch_fit_to_data_frame <- function(model, data) {
  records <- vector("list", length(data$catch$frames))
  start <- 1L
  for (series in seq_along(data$catch$frames)) {
    frame <- data$catch$frames[[series]]
    rows <- nrow(frame)
    idx <- start:(start + rows - 1L)
    records[[series]] <- data.frame(
      series = series,
      year = as.integer(frame$year),
      fleet = data$fleet_names[as.integer(frame$fleet)],
      season = as.integer(frame$season),
      sex = c("both", "male", "female")[as.integer(frame$sex) + 1L],
      predicted = as.numeric(model$pre_catch[idx]),
      residual = as.numeric(model$res_catch[idx]),
      stringsAsFactors = FALSE
    )
    start <- start + rows
  }
  do.call(rbind, records)
}

compare_catch_fit_to_admb <- function(model, inputs) {
  data <- make_gmacs_rtmb_data(inputs)
  rtmb <- catch_fit_to_data_frame(model, data)
  admb <- inputs$admb_catch_fit
  merged <- merge(
    admb,
    rtmb,
    by = c("series", "year", "fleet", "season", "sex"),
    suffixes = c("_admb", "_rtmb")
  )
  merged$pred_abs_diff <- abs(merged$predicted_admb - merged$predicted_rtmb)
  merged$resid_abs_diff <- abs(merged$residual_admb - merged$residual_rtmb)
  list(
    predicted_max_abs_diff = max(merged$pred_abs_diff),
    residual_max_abs_diff = max(merged$resid_abs_diff),
    n = nrow(merged),
    data = merged
  )
}

index_fit_to_data_frame <- function(model, data) {
  survey <- data$survey$data
  data.frame(
    series = as.integer(survey$index),
    year = as.integer(survey$year),
    fleet = data$fleet_names[as.integer(survey$fleet)],
    season = as.integer(survey$season),
    sex = c("both", "male", "female")[as.integer(survey$sex) + 1L],
    maturity = c("all", "mature", "immature")[as.integer(survey$maturity) + 1L],
    predicted = as.numeric(model$pre_cpue),
    residual = as.numeric(model$res_cpue),
    q = as.numeric(model$survey_q_calc[as.integer(survey$index)]),
    stringsAsFactors = FALSE
  )
}

compare_index_fit_to_admb <- function(model, inputs) {
  data <- make_gmacs_rtmb_data(inputs)
  rtmb <- index_fit_to_data_frame(model, data)
  admb <- inputs$admb_index_fit
  merged <- merge(
    admb,
    rtmb,
    by = c("series", "year", "fleet", "season", "sex", "maturity"),
    suffixes = c("_admb", "_rtmb")
  )
  merged$pred_abs_diff <- abs(merged$predicted_admb - merged$predicted_rtmb)
  merged$q_abs_diff <- abs(merged$q_admb - merged$q_rtmb)
  list(
    predicted_max_abs_diff = max(merged$pred_abs_diff),
    q_max_abs_diff = max(merged$q_abs_diff),
    n = nrow(merged),
    data = merged
  )
}

size_composition_to_data_frame <- function(model, data) {
  controls <- data$size_comp_controls
  records <- vector("list", max(controls$comp_aggregator))
  start <- 1L
  for (modified in seq_len(max(controls$comp_aggregator))) {
    series <- which(controls$comp_aggregator == modified)[1]
    frame <- data$size_comp$frames[[series]]
    rows <- model$size_comp_rows[modified]
    cols <- model$size_comp_cols[modified]
    idx <- start:(start + rows * cols - 1L)
    pred <- matrix(as.numeric(model$pre_size_comps[idx]), nrow = rows, ncol = cols, byrow = TRUE)
    records[[modified]] <- data.frame(
      modified_series = modified,
      year = rep(as.integer(frame$year), each = ncol(pred)),
      fleet = rep(data$fleet_names[as.integer(frame$fleet)], each = ncol(pred)),
      season = rep(as.integer(frame$season), each = ncol(pred)),
      size_bin = rep(seq_len(ncol(pred)), times = nrow(pred)),
      value = as.vector(t(pred)),
      stringsAsFactors = FALSE
    )
    start <- start + rows * cols
  }
  do.call(rbind, records)
}

compare_size_composition_to_admb <- function(model, inputs) {
  data <- make_gmacs_rtmb_data(inputs)
  rtmb <- size_composition_to_data_frame(model, data)
  admb <- inputs$admb_size_fit[inputs$admb_size_fit$vector == "pred", ]
  merged <- merge(
    admb,
    rtmb,
    by = c("modified_series", "year", "fleet", "season", "size_bin"),
    suffixes = c("_admb", "_rtmb")
  )
  merged$abs_diff <- abs(merged$value_admb - merged$value_rtmb)
  list(
    max_abs_diff = max(merged$abs_diff),
    mean_abs_diff = mean(merged$abs_diff),
    n = nrow(merged),
    data = merged
  )
}

validate_gmacs_inputs <- function(inputs) {
  dat <- inputs$data
  pars <- inputs$parameters
  dims <- dat$dimensions
  checks <- c(
    nclass_size_breaks = length(dat$size_breaks) == dims[["nclass"]] + 1L,
    n_size_sex = length(dat$n_size_sex) == dims[["nsex"]],
    season_indices = all(c(dat$season_recruitment, dat$season_growth, dat$season_N) >= 1L) &&
      all(c(dat$season_recruitment, dat$season_growth, dat$season_N) <= dims[["nseason"]]),
    season_type = length(dat$season_type) == dims[["nseason"]],
    catch_rows = sum(dat$catch$rows) == sum(vapply(dat$catch$frames, nrow, integer(1))),
    survey_rows = nrow(dat$survey$data) > 0L,
    size_comp_rows = sum(dat$size_comp$rows) ==
      sum(vapply(dat$size_comp$frames, nrow, integer(1))),
    theta_present = "theta" %in% names(pars),
    growth_present = "Grwth" %in% names(pars),
    selectivity_present = "log_slx_pars" %in% names(pars)
  )
  if (!all(checks)) {
    stop(
      "Input validation failed: ",
      paste(names(checks)[!checks], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(checks)
}

gmacs_rtmb_nll_factory <- function(data) {
  force(data)

  function(par) {
    nll <- 0
    model <- initialize_bbrkc_model_parameters(par, data, predictions = TRUE)

    ## Keep every parameter block on the AD tape while deterministic pieces are
    ## ported one at a time. This is intentionally neutral in the objective.
    for (nm in names(par)) {
      nll <- nll + sum(par[[nm]] * 0)
    }

    zero_ad <- sum(par$theta * 0)
    nlog_penalty <- calc_bbrkc_penalties(par, data, model, zero_ad)
    weighted_nlog_penalty <- nlog_penalty * data$likelihood_emphasis$penalty
    prior_density <- calc_bbrkc_prior_density(par, data, zero_ad)

    catch_likelihood <- function() {
      out <- vector("list", length(data$catch$frames))
      start <- 1L
      for (series in seq_along(data$catch$frames)) {
        frame <- data$catch$frames[[series]]
        rows <- nrow(frame)
        idx <- start:(start + rows - 1L)
        catch_sd <- sqrt(log(1 + frame$cv^2))
        out[[series]] <- admb_dnorm_nll(model$res_catch[idx], catch_sd)
        start <- start + rows
      }
      do.call(c, out)
    }

    index_likelihood <- function() {
      survey <- data$survey$data
      controls <- data$q_controls
      out <- vector("list", length(data$survey$index_type))
      for (series in seq_along(out)) {
        rows <- which(survey$index == series)
        nll_series <- zero_ad
        for (row in rows) {
          cvadd2 <- if (controls$add_cv_links[series] > 0L) {
            log(1 + exp(par$log_add_cv[controls$add_cv_links[series]])^2)
          } else {
            0
          }
          cvobs2 <- log(1 + survey$cv[row]^2) / controls$cpue_lambda[series]
          stdtmp <- sqrt(cvobs2 + cvadd2)
          nll_series <- nll_series + log(stdtmp) +
            0.5 * (model$res_cpue[row] / stdtmp)^2
        }
        out[[series]] <- nll_series
      }
      do.call(c, out)
    }

    length_likelihood <- function() {
      obs <- data$observed_size_compositions
      out <- vector("list", length(obs$observed))
      start <- 1L
      for (series in seq_along(obs$observed)) {
        rows <- nrow(obs$observed[[series]])
        cols <- ncol(obs$observed[[series]])
        idx <- start:(start + rows * cols - 1L)
        pred <- model$pre_size_comps[idx]
        log_effn <- par$log_vn[series] + log(obs$sample_size[[series]] * obs$lambda[series])
        out[[series]] <- switch(
          as.character(obs$likelihood_type[series]),
          "0" = zero_ad,
          "1" = multinomial_nll_flat(obs$observed[[series]], pred, log_effn),
          "2" = robust_multinomial_nll_flat(obs$observed[[series]], pred, log_effn),
          "5" = dirichlet_nll_flat(obs$observed[[series]], pred, log_effn),
          stop("Unsupported size-composition likelihood type: ", obs$likelihood_type[series], call. = FALSE)
        )
        start <- start + rows * cols
      }
      do.call(c, out)
    }

    recruitment_likelihood <- function() {
      sigR <- exp(par$theta[10])
      out <- vector("list", 3)
      out[[1]] <- admb_dnorm_nll(model$res_recruit, sigR)
      if (length(par$rec_ini) > 0L && data$rec_controls$rec_ini_phz > 0) {
        out[[2]] <- admb_dnorm_nll(par$rec_ini, sigR)
      } else {
        out[[2]] <- zero_ad
      }
      if (length(par$logit_rec_prop_est) > 0L && data$rec_controls$rec_prop_phz > 0) {
        out[[3]] <- admb_dnorm_nll(par$logit_rec_prop_est, 2)
      } else {
        out[[3]] <- zero_ad
      }
      do.call(c, out)
    }

    growth_likelihood <- function() {
      ## BBRKC has no active growth observations in the current data file.
      zero_ad
    }

    catch_nloglike <- catch_likelihood()
    index_nloglike <- index_likelihood()
    length_nloglike <- length_likelihood()
    recruitment_nloglike <- recruitment_likelihood()
    growth_nloglike <- growth_likelihood()

    nloglike <- c(
      catch = sum(catch_nloglike * data$likelihood_emphasis$catch),
      index = sum(index_nloglike * data$likelihood_emphasis$cpue_emphasis),
      length = sum(length_nloglike * data$observed_size_compositions$emphasis),
      recruitment = sum(recruitment_nloglike),
      growth = growth_nloglike
    )

    RTMB::REPORT(nloglike)
    RTMB::REPORT(catch_nloglike)
    RTMB::REPORT(index_nloglike)
    RTMB::REPORT(length_nloglike)
    RTMB::REPORT(recruitment_nloglike)
    RTMB::REPORT(growth_nloglike)
    RTMB::REPORT(nlog_penalty)
    RTMB::REPORT(weighted_nlog_penalty)
    RTMB::REPORT(prior_density)
    RTMB::REPORT(data$dimensions)
    RTMB::REPORT(data$mid_points)
    RTMB::REPORT(model$M0)
    RTMB::REPORT(model$logRbar)
    RTMB::REPORT(model$log_slx_capture)
    RTMB::REPORT(model$log_slx_retaind)
    RTMB::REPORT(model$log_slx_discard)
    RTMB::REPORT(model$M)
    RTMB::REPORT(model$ft)
    RTMB::REPORT(model$fout)
    RTMB::REPORT(model$F)
    RTMB::REPORT(model$F2)
    RTMB::REPORT(model$Z)
    RTMB::REPORT(model$Z2)
    RTMB::REPORT(model$S)
    RTMB::REPORT(model$molt_increment)
    RTMB::REPORT(model$gscale)
    RTMB::REPORT(model$molt_probability)
    RTMB::REPORT(model$growth_transition)
    RTMB::REPORT(model$rec_sdd)
    RTMB::REPORT(model$recruits)
    RTMB::REPORT(model$rec_dev)
    RTMB::REPORT(model$res_recruit)
    RTMB::REPORT(model$logit_rec_prop)
    RTMB::REPORT(model$logN0)
    RTMB::REPORT(model$N0)
    RTMB::REPORT(model$N)
    RTMB::REPORT(model$pre_catch)
    RTMB::REPORT(model$res_catch)
    RTMB::REPORT(model$log_q_catch)
    RTMB::REPORT(model$vulnerable_cpue)
    RTMB::REPORT(model$pre_cpue)
    RTMB::REPORT(model$res_cpue)
    RTMB::REPORT(model$survey_q_calc)
    RTMB::REPORT(model$pre_size_comps)
    RTMB::REPORT(model$size_comp_rows)
    RTMB::REPORT(model$size_comp_cols)

    nll + sum(nloglike) + sum(weighted_nlog_penalty) + sum(prior_density)
  }
}

make_gmacs_rtmb_object <- function(inputs, random = NULL, map = make_bbrkc_rtmb_map(inputs)) {
  if (!requireNamespace("RTMB", quietly = TRUE)) {
    stop("The RTMB package is required to tape the objective.", call. = FALSE)
  }
  validate_gmacs_inputs(inputs)
  data <- make_gmacs_rtmb_data(inputs)
  parameters <- make_gmacs_parameter_list(inputs)
  RTMB::MakeADFun(gmacs_rtmb_nll_factory(data), parameters, random = random, map = map)
}

# GMACS-R

This repository contains an R/RTMB scaffold for porting GMACS stock assessment
logic, with the current implementation focused on the BBRKC model.

The detailed porting notes, parity checks, and next implementation steps live in
`bbrkc_rtmb_port_steps.qmd`.

## Repository contents

- `R/read_gmacs_bbrkc.R`: parsers for BBRKC GMACS input files and ADMB reference
  outputs.
- `R/gmacs_rtmb_nll.R`: RTMB data/parameter construction, parameter mapping,
  deterministic model blocks, comparison helpers, and the RTMB objective
  factory.
- `run_bbrkc.R`: parser, RTMB tape, and component-comparison smoke run.
- `optimize_bbrkc.R`: short `nlminb()` optimization smoke test.
- `run_tanner.R`: Tanner input parser and fitted-parameter summary for the
  files in `examples/tanners`; it also discovers Tanner ADMB report outputs
  when `Gmacsall.out` and `gmacs.rep` or recognized equivalent names are added.
- `run_snow.R`: Snow crab input parser and ADMB report-summary smoke run for
  the files in `examples/snow`.
- `gmacs.tpl`: GMACS ADMB template reference.
- `examples/tanners/`: Tanner crab GMACS input files retained as local
  reference material.
- `GMACS_InputFiles.zip`: original archive for the Tanner crab input files. This
  is ignored for future commits; the extracted `tanners/` files are the useful
  repository copy.

## Inputs

The BBRKC scripts require these files in one directory:

- `gmacs.dat`
- `gmacs.pin`
- `Gmacsall.out`
- `gmacs.rep`

By default the scripts look in `examples/BBRKC`, which contains the committed
BBRKC parity example. To use another location, set `GMACS_BBRKC_ROOT`:

```sh
GMACS_BBRKC_ROOT=/path/to/BBRKC Rscript run_bbrkc.R
GMACS_BBRKC_ROOT=/path/to/BBRKC Rscript optimize_bbrkc.R
```

From R:

```r
Sys.setenv(GMACS_BBRKC_ROOT = "/path/to/BBRKC")
source("run_bbrkc.R")
```

## Current scope

- Reads BBRKC data sections for dimensions, M proportions, fleets, catch,
  survey indices, size compositions, growth, and environment.
- Reads fitted parameter vectors from `gmacs.pin`.
- Builds an RTMB `map` from BBRKC control phases, mapping the 560 full `.pin`
  entries down to the 366 active ADMB parameters.
- Reads ADMB likelihood, selectivity, fishing mortality, natural mortality,
  total mortality, growth, numbers-at-size, catch-fit, index-fit, size-fit, and
  summary reference blocks from `Gmacsall.out` and `gmacs.rep`.
- Ports deterministic BBRKC selectivity, fishing mortality, natural mortality,
  total mortality, survival setup, molt probability, growth transition matrices,
  recruitment, FREEPARSSCALED initial numbers, seasonal population updates,
  catch predictions, survey index predictions, and size-composition
  predictions.
- Computes gamma growth transitions and recruitment size distributions with
  `RTMB::pgamma()` so active growth and recruitment-distribution parameters
  remain on the AD tape.
- Builds an RTMB objective with named likelihood slots matching `gmacs.tpl`.
- Audits the active likelihood components, raw `nlogPenalty`, and
  `priorDensity` against ADMB reference blocks in `run_bbrkc.R`.

## Rendered Notes

The Quarto notes are rendered for GitHub Pages under `docs/`:

```sh
quarto render bbrkc_rtmb_port_steps.qmd --output-dir docs --output bbrkc.html --no-clean
quarto render tanner_rtmb_port_steps.qmd --output-dir docs --output tanner.html --no-clean
quarto render snow_rtmb_port_steps.qmd --output-dir docs --output snow.html --no-clean
```

Configure GitHub Pages to serve the `docs/` directory from the `main` branch.

## Next Porting Steps

1. Run `Rscript run_bbrkc.R` after every deterministic prediction change and
   inspect the printed ADMB likelihood/penalty/prior audit.
2. Add active growth likelihood terms if future GMACS inputs include growth
   observations.
3. Broaden validation against additional GMACS model configurations.

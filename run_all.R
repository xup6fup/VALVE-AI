# VALVE-AI: generate 21 PDF figures and five CSV tables from the prepared data.
# RStudio: open VALVE-AI.Rproj, then run source("run_all.R").
# Command line: Rscript run_all.R [input_csv] [output_directory]
.valve_project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(.valve_project_root, "R", "01_data.R")))
  stop("Run this script from the VALVE-AI project directory.")
.valve_args <- commandArgs(trailingOnly = TRUE)
.valve_input <- if (length(.valve_args)) .valve_args[1] else
  Sys.getenv("VALVE_AI_DATA", unset = file.path(.valve_project_root, "data", "VALVE_AI.csv"))
.valve_output <- if (length(.valve_args) > 1L) .valve_args[2] else
  file.path(.valve_project_root, "output")
if (!file.exists(.valve_input))
  stop("Place the approved analysis dataset at data/VALVE_AI.csv or specify its path. See README.md.")
for (module in c("01_data.R", "02_statistics.R", "03_tables.R", "04_figures.R",
                 "05_exports.R", "06_tte_sensitivity.R"))
  source(file.path(.valve_project_root, "R", module), encoding = "UTF-8")
result <- run_valve_ai(.valve_input, .valve_output)

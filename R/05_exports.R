
valve_draw_figure <- function(id, df, config) {
  if (id == "Ex_Data_Fig_01") valve_draw_ecg(df, config)
  else if (id %in% c(sprintf("Fig_%02d", 2:5), sprintf("Ex_Data_Fig_%02d", 2:4)))
    valve_draw_forest(df, id, config)
  else valve_draw_yield(df, id, config)
}

valve_dimensions <- function(id, df) {
  if (id == "Ex_Data_Fig_01") c(width = 183 / 25.4, height = 124.037 / 25.4)
  else if (id %in% c(sprintf("Fig_%02d", 2:5), sprintf("Ex_Data_Fig_%02d", 2:4)))
    valve_forest_dimensions(df, id)
  else valve_yield_dimensions(df, id)
}

valve_save_figure <- function(id, df, path, config) {
  size <- valve_dimensions(id, df)
  grDevices::cairo_pdf(path, width = size["width"], height = size["height"],
                       family = config$font_family, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  valve_draw_figure(id, df, config)
  invisible(path)
}

run_valve_ai <- function(input, output_dir = "output", config = valve_config()) {
  if (!capabilities("cairo")) stop("This R installation requires Cairo PDF support.")
  if (!requireNamespace("logistf", quietly = TRUE))
    stop("Install the analysis dependency first: install.packages('logistf')")

  message("[1/4] Read the prepared analysis dataset")
  d <- valve_read_input(input, config)
  message("[2/4] Calculate figure estimates and the five tables")
  figures <- c(list(Ex_Data_Fig_01 = valve_ecg_source(d)), valve_analyse(d, config))
  ids <- c(sprintf("Fig_%02d", 2:5), sprintf("Ex_Data_Fig_%02d", 1:10),
           sprintf("Supplemental_Fig_%02d", 1:8))
  tables <- valve_tables(d, config)

  figure_dir <- file.path(output_dir, "figures")
  table_dir <- file.path(output_dir, "tables")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  message("[3/4] Write five CSV tables")
  for (id in names(tables))
    write.csv(tables[[id]], file.path(table_dir, paste0(id, ".csv")),
              row.names = FALSE, na = "", fileEncoding = "UTF-8")

  message("[4/4] Draw 22 individual PDF figures")
  for (id in ids) {
    message("  ", id)
    valve_save_figure(id, figures[[id]], file.path(figure_dir, paste0(id, ".pdf")), config)
  }
  message("Done: 22 PDF figures and 5 CSV tables in ", output_dir)
  invisible(output_dir)
}

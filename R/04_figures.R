
valve_sg_spec <- function(key) {
  spec <- VALVE_SUBGROUPS[[key]]
  if (is.null(spec)) stop("Missing subgroup metadata: ", key)
  if (is.null(spec$display_levels)) spec$display_levels <- spec$levels
  if (length(spec$display_levels) != length(spec$levels))
    stop("Subgroup levels/display_levels differ: ", key)
  spec
}

valve_sg_values <- function(d, key) {
  spec <- valve_sg_spec(key)
  field <- if (is.null(spec$column)) key else spec$column
  if (!field %in% names(d)) stop("Missing subgroup input: ", field)
  values <- as.character(d[[field]])
  unexpected <- unique(values[!is.na(values) & !values %in% spec$levels])
  if (length(unexpected)) stop("Unrecognized subgroup values in ", field, ": ", paste(unexpected, collapse = " | "))
  values
}

valve_outcome <- function(key) {
  spec <- VALVE_OUTCOMES[[key]]
  if (is.null(spec)) stop("Missing outcome metadata: ", key)
  spec
}

valve_proportions <- function(d, figure, variables = character(),
                              outcome = "outcome_vhd", conditional = FALSE,
                              include_all = FALSE, config = list()) {
  # Conditional analyses restrict the denominator to participants with detected
  # VHD; all other yield analyses use TTE recipients. Each outcome still uses
  # its own non-missing denominator (e.g. cardiology-referral eligibility).
  base <- d[!is.na(d$PATHWAY_INDEX), , drop = FALSE]
  if (conditional) base <- base[base$outcome_vhd %in% 1, , drop = FALSE]
  specs <- lapply(variables, valve_sg_spec); names(specs) <- variables
  if (include_all) specs <- c(list(All = list(name = "All participants",
                levels = "All participants", display_levels = "All participants")), specs)
  rows <- list(); outcome_meta <- valve_outcome(outcome)
  for (variable in names(specs)) {
    spec <- specs[[variable]]
    values <- if (variable == "All") rep("All participants", nrow(base)) else valve_sg_values(base, variable)
    for (j in seq_along(spec$levels)) {
      selected <- base[values %in% spec$levels[j], , drop = FALSE]
      cells <- lapply(1:5, function(path) selected[selected$PATHWAY_INDEX %in% path, , drop = FALSE])
      counts <- lapply(cells, function(z) c(Events = sum(z[[outcome]] == 1, na.rm = TRUE),
                                          N = sum(!is.na(z[[outcome]]))))
      re <- counts[[1]]["Events"]; rn <- counts[[1]]["N"]
      for (path in 1:5) {
        e <- unname(counts[[path]]["Events"]); n <- unname(counts[[path]]["N"])
        interval <- valve_wilson(e, n)
        test <- if (path == 1L) list(P = NA_real_, Method = "Reference", Status = if (n) "Reference" else "NE: empty reference denominator", Warnings = "") else
          valve_fisher(rbind(c(e, n-e), c(re, rn-re)), config)
        rows[[length(rows)+1L]] <- data.frame(
          Figure = figure, Variable = variable, Subgroup = spec$name,
          Level = spec$display_levels[j], Outcome = outcome_meta$label,
          Outcome_variable = outcome, Pathway_index = path, Pathway = VALVE_PATHWAYS[path],
          Eligible_records = nrow(cells[[path]]), Events = e, N = n,
          Missing_outcome = nrow(cells[[path]])-n,
          Missing_subgroup = if (variable == "All") 0L else sum(is.na(values) | !values %in% spec$levels),
          Proportion = interval["Proportion"], Lower = interval["Lower"], Upper = interval["Upper"],
          P = test$P, P_display = if (path == 1L && n > 0) "Ref." else valve_p_display(test$P),
          Reference = VALVE_PATHWAYS[1], Reference_events = unname(re), Reference_N = unname(rn),
          CI_method = "Two-sided Wilson score 95% CI; no continuity correction",
          Test_method = "Two-sided Fisher exact test", Test_status = test$Status,
          Test_warnings = test$Warnings,
          Ascertainment = "Supplied 90-day indicators; recorded TTE pathway",
          stringsAsFactors = FALSE, row.names = NULL)
      }
    }
  }
  valve_rbind(rows)
}

valve_analyse <- function(d, config = list()) {
  clinical <- c("subgroup_sex", "subgroup_age", "subgroup_recent_relevant_symptoms", "subgroup_outpatient_cardiology_care", "subgroup_cardiovascular_prescription", "subgroup_diabetes_mellitus_hypertension_histories", "subgroup_chronic_kidney_disease_history", "subgroup_cardiovascular_disease_history")
  ecg <- c("subgroup_most_recent_ecg_ai_risk", "subgroup_time_since_most_recent_ai_high_risk_ecg", "subgroup_ai_ecg_risk_profile_preceding_year")
  condition <- "prespecified_clinical_triggers_for_tte"
  # These keys identify the supplied variables, not newly reconstructed clinical
  # definitions. The plotting metadata controls display labels and category order.
  d$randomization_group <- factor(d$randomization_group, levels = c("Control", "Intervention"))
  high <- d[d$AI_HIGH %in% TRUE, , drop = FALSE]
  primary_secondary <- c("outcome_vhd", c("outcome_aortic_stenosis", "outcome_aortic_regurgitation", "outcome_mitral_regurgitation", "outcome_tricuspid_regurgitation", "outcome_other_structural_abnormality"))
  other <- c("outcome_vhd_with_relevant_symptoms", "outcome_vhd_without_relevant_symptoms", "outcome_vhd_ntprobnp_ge900", "outcome_left_sided_vhd", "outcome_left_sided_vhd_lvedd_ge60", "outcome_severe_vhd", "outcome_severe_vhd_surgical_criteria")
  exam <- c("outcome_vhd_composite_examinations", "outcome_vhd_chest_xray", "outcome_vhd_ntprobnp_test", "outcome_vhd_high_sensitivity_troponin_i_test", "outcome_vhd_exercise_test", "outcome_vhd_transesophageal_echo")
  management <- c("outcome_vhd_composite_managements", "outcome_vhd_cardiologist_referral", "outcome_vhd_pharmacotherapy_modification", "outcome_vhd_hf_pharmacotherapy_modification", "outcome_vhd_antihypertensive_pharmacotherapy_modification", "outcome_vhd_cardiovascular_surgeon_referral", "outcome_vhd_valve_intervention")
  requests <- list(); interactions <- list()

  add_effect <- function(figure, panel, population, q, outcome, block,
                         section = NULL, variable = "", level = "", interact = FALSE) {
    meta <- valve_outcome(outcome); counts <- valve_counts(q, outcome)
    key <- paste(figure, panel, outcome, variable, level, sep = ":")
    req <- data.frame(Figure = figure, Panel = panel, Population = population,
      Section = if (is.null(section)) meta$section else section,
      Outcome = meta$label, Outcome_variable = outcome,
      Subgroup = if (variable %in% names(VALVE_SUBGROUPS)) valve_sg_spec(variable)$name else if (variable == "AI_RISK") "AI-ECG risk" else variable,
      Variable = variable, Level = level, Request = key, Block = block,
      N_intervention = counts$N[1], Events_intervention = counts$Events[1],
      N_control = counts$N[2], Events_control = counts$Events[2],
      Missing_intervention = counts$Missing[1], Missing_control = counts$Missing[2],
      Available = TRUE, Unavailable_reason = "", Proxy = FALSE, Proxy_reason = "",
      Interaction_request = if (interact) block else "", stringsAsFactors = FALSE)
    requests[[length(requests)+1L]] <<- list(row = req, counts = counts)
  }
  add_subgroup <- function(figure, population, q, variable, outcomes) {
    spec <- valve_sg_spec(variable); values <- valve_sg_values(q, variable)
    for (outcome in outcomes) {
      block <- paste(figure, variable, outcome, sep = ":")
      cnt <- list()
      for (j in seq_along(spec$levels)) {
        cell <- q[values %in% spec$levels[j], , drop = FALSE]
        add_effect(figure, "", population, cell, outcome, block, section = spec$name,
                   variable = variable, level = spec$display_levels[j], interact = TRUE)
        z <- valve_counts(cell, outcome); z$Stratum <- spec$display_levels[j]; cnt[[j]] <- z
      }
      interactions[[block]] <<- do.call(rbind, cnt)
    }
  }

  for (outcome in c(primary_secondary, "outcome_new_psychiatric_consultation"))
    add_effect("Fig_02", "", "All participants", d, outcome, paste0("Fig_02:", outcome))
  for (outcome in c("outcome_tte_after_randomization", primary_secondary)) {
    block <- paste0("Fig_03:", outcome); cnt <- list()
    for (risk in c("High", "Low")) {
      q <- d[d$AI_HIGH %in% (risk == "High"), , drop = FALSE]
      add_effect("Fig_03", "", paste("AI-ECG", tolower(risk), "risk"), q,
        outcome, block, section = valve_outcome(outcome)$label,
        variable = "AI_RISK", level = risk, interact = TRUE)
      z <- valve_counts(q, outcome); z$Stratum <- risk; cnt[[risk]] <- z
    }
    interactions[[block]] <- do.call(rbind, cnt)
  }
  for (figure in c("Fig_04", "Fig_05")) {
    outcomes <- if (figure == "Fig_04") other else c(exam, management)
    for (outcome in outcomes) {
      block <- paste(figure, outcome, sep = ":")
      add_effect(figure, "a", "All participants", d, outcome, block)
      add_effect(figure, "b", "AI-ECG high risk", high, outcome, block)
    }
  }
  for (variable in c(condition, clinical)) {
    add_subgroup("Ex_Data_Fig_02", "All participants", d, variable, "outcome_vhd")
    add_subgroup("Ex_Data_Fig_03", "AI-ECG high risk", high, variable, "outcome_vhd")
  }
  for (variable in ecg)
    add_subgroup("Ex_Data_Fig_04", "AI-ECG high risk", high, variable, c("outcome_tte_after_randomization", "outcome_vhd"))

  # First determine the fitting method for each complete, coherent block.
  # A zero-event contrast affects every estimable effect and interaction within
  # that outcome/population/subgroup block. Both-arm-zero contrasts stay NE.
  blocks <- unique(vapply(requests, function(z) z$row$Block, character(1)))
  modes <- setNames(rep("glm", length(blocks)), blocks)
  for (req in requests) {
    z <- req$counts
    if (all(z$N > 0) && sum(z$Events == 0) == 1L) modes[req$row$Block] <- "Firth"
  }
  interaction_results <- lapply(names(interactions), function(block)
    valve_interaction(interactions[[block]], modes[[block]]))
  names(interaction_results) <- names(interactions)
  rows <- lapply(requests, function(req) {
    result <- valve_effect(req$counts, modes[[req$row$Block]])
    ix <- interaction_results[[req$row$Interaction_request]]
    if (is.null(ix)) ix <- list(Method = NA_character_, Status = NA_character_, P = NA_real_,
                               LR = NA_real_, DF = NA_real_, Warnings = NA_character_)
    names(ix) <- paste0("Interaction_", names(ix))
    row <- cbind(req$row, as.data.frame(result, stringsAsFactors = FALSE),
                 as.data.frame(ix, stringsAsFactors = FALSE))
    row$P_display <- valve_p_display(row$P, row$Method)
    row$Interaction_P_display <- valve_p_display(row$Interaction_P, row$Interaction_Method)
    row
  })
  effects <- valve_rbind(rows)
  # ECG forest follows variable -> category -> TTE, primary outcome.
  ix <- which(effects$Figure == "Ex_Data_Fig_04")
  if (length(ix)) {
    z <- effects[ix, , drop = FALSE]
    lev_rank <- mapply(function(v, lev) match(lev, valve_sg_spec(v)$display_levels), z$Variable, z$Level)
    z <- z[order(match(z$Variable, ecg), lev_rank, match(z$Outcome_variable, c("outcome_tte_after_randomization", "outcome_vhd"))), , drop = FALSE]
    effects[ix, ] <- z
  }
  output <- split(effects, factor(effects$Figure, levels = unique(effects$Figure)))
  output <- lapply(output, function(z) { rownames(z) <- NULL; z })

  c(output, valve_analyse_yields(d, config))
}

valve_analyse_yields <- function(d, config = list()) {
  clinical <- c("subgroup_sex", "subgroup_age", "subgroup_recent_relevant_symptoms", "subgroup_outpatient_cardiology_care", "subgroup_cardiovascular_prescription", "subgroup_diabetes_mellitus_hypertension_histories", "subgroup_chronic_kidney_disease_history", "subgroup_cardiovascular_disease_history")
  ecg <- c("subgroup_most_recent_ecg_ai_risk", "subgroup_time_since_most_recent_ai_high_risk_ecg", "subgroup_ai_ecg_risk_profile_preceding_year")
  condition <- "prespecified_clinical_triggers_for_tte"
  other <- c("outcome_vhd_with_relevant_symptoms", "outcome_vhd_without_relevant_symptoms", "outcome_vhd_ntprobnp_ge900", "outcome_left_sided_vhd", "outcome_left_sided_vhd_lvedd_ge60", "outcome_severe_vhd", "outcome_severe_vhd_surgical_criteria")
  exam <- c("outcome_vhd_composite_examinations", "outcome_vhd_chest_xray", "outcome_vhd_ntprobnp_test", "outcome_vhd_high_sensitivity_troponin_i_test", "outcome_vhd_exercise_test", "outcome_vhd_transesophageal_echo")
  management <- c("outcome_vhd_composite_managements", "outcome_vhd_cardiologist_referral", "outcome_vhd_pharmacotherapy_modification", "outcome_vhd_hf_pharmacotherapy_modification", "outcome_vhd_antihypertensive_pharmacotherapy_modification", "outcome_vhd_cardiovascular_surgeon_referral", "outcome_vhd_valve_intervention")
  output <- list()

  # Disease-detection tree (Extended Data Fig. 5). Unlike yield comparisons,
  # layer 1 and layer 2 denominators include all randomized participants.
  tree <- list()
  tree_row <- function(q, arm, risk, source, layer) {
    e <- sum(q$outcome_vhd == 1, na.rm = TRUE); n <- sum(!is.na(q$outcome_vhd))
    w <- valve_wilson(e, n)
    data.frame(randomization_group = arm, AI_RISK = risk, SOURCE = source, Layer = layer,
      Events = e, Total = nrow(q), Observed = n,
      CI_lower = w["Lower"], CI_upper = w["Upper"], stringsAsFactors = FALSE, row.names = NULL)
  }
  for (arm in c("Intervention", "Control")) {
    a <- d[d$randomization_group %in% arm, , drop = FALSE]
    tree[[length(tree)+1L]] <- tree_row(a, arm, "All", "All", 1L)
    for (risk in c("High", "Low")) {
      q <- a[a$AI_HIGH %in% (risk == "High"), , drop = FALSE]
      tree[[length(tree)+1L]] <- tree_row(q, arm, risk, "All", 2L)
      sources <- if (arm == "Intervention" && risk == "High") c("Trial-triggered", "Physician-triggered", "No TTE") else c("Physician-triggered", "No TTE")
      for (source in sources) {
        keep <- if (source == "Trial-triggered") q$PATHWAY_INDEX %in% 1 else
          if (source == "Physician-triggered") q$PATHWAY_INDEX %in% 2:5 else q$outcome_tte_after_randomization %in% 0
        tree[[length(tree)+1L]] <- tree_row(q[keep, , drop = FALSE], arm, risk, source, 3L)
      }
      unknown <- is.na(q$outcome_tte_after_randomization) | (q$outcome_tte_after_randomization %in% 1 & is.na(q$PATHWAY_INDEX))
      if (any(unknown)) tree[[length(tree)+1L]] <- tree_row(q[unknown, , drop = FALSE], arm, risk, "Unknown TTE status/pathway", 3L)
    }
    if (anyNA(a$AI_HIGH)) tree[[length(tree)+1L]] <- tree_row(a[is.na(a$AI_HIGH), , drop = FALSE], arm, "Unknown", "All", 2L)
  }
  tree_rows <- valve_rbind(tree)
  tree_rows$Figure <- "Ex_Data_Fig_05"; tree_rows$Panel <- "a"
  overall_yield <- valve_proportions(d, "Ex_Data_Fig_05", include_all = TRUE, config = config)
  overall_yield$Panel <- "b"
  output$Ex_Data_Fig_05 <- valve_rbind(list(tree_rows, overall_yield))

  # Extended Data Fig. 6: one overall exact yield test across subgroup levels,
  # exclusively within trial-triggered TTE. Multilevel variables get one P.
  trial <- d[d$PATHWAY_INDEX %in% 1, , drop = FALSE]; trial_rows <- list()
  variables <- c(ecg, condition, clinical)
  for (vi in seq_along(variables)) {
    variable <- variables[vi]; spec <- valve_sg_spec(variable); values <- valve_sg_values(trial, variable)
    cells <- lapply(spec$levels, function(lev) trial[values %in% lev, , drop = FALSE])
    es <- vapply(cells, function(z) sum(z$outcome_vhd == 1, na.rm = TRUE), integer(1))
    ns <- vapply(cells, function(z) sum(!is.na(z$outcome_vhd)), integer(1))
    test <- valve_fisher(cbind(es, ns-es), config)
    for (j in seq_along(cells)) {
      w <- valve_wilson(es[j], ns[j])
      trial_rows[[length(trial_rows)+1L]] <- data.frame(Figure = "Ex_Data_Fig_06",
        Test_ID = paste0("trial_triggered:", variable), Variable = variable, Subgroup = spec$name,
        Variable_order = vi, Level = spec$display_levels[j], Level_order = j,
        Pathway_index = 1L, Pathway = VALVE_PATHWAYS[1], Outcome = valve_outcome("outcome_vhd")$label,
        Outcome_variable = "outcome_vhd", Events = es[j], N = ns[j], Non_events = ns[j]-es[j],
        Total_TTE = nrow(cells[[j]]), Missing_outcome = nrow(cells[[j]])-ns[j],
        Missing_subgroup = sum(is.na(values) | !values %in% spec$levels),
        Proportion = w["Proportion"], Lower = w["Lower"], Upper = w["Upper"],
        n_N = paste0(es[j], "/", ns[j]),
        Proportion_display = if (ns[j]) sprintf("%.1f%%", 100*w["Proportion"]) else "NE",
        CI_display = if (ns[j]) sprintf("%.1f%%-%.1f%%", 100*w["Lower"], 100*w["Upper"]) else "NE",
        CI_method = "Two-sided Wilson score 95% CI without continuity correction",
        P = test$P, Pglobal = test$P, P_display = valve_p_display(test$P),
        Test_method = test$Method, Test_status = test$Status, Test_warnings = test$Warnings,
        stringsAsFactors = FALSE, row.names = NULL)
    }
  }
  output$Ex_Data_Fig_06 <- valve_rbind(trial_rows)

  secondary <- c("outcome_aortic_stenosis", "outcome_aortic_regurgitation", "outcome_mitral_regurgitation", "outcome_tricuspid_regurgitation", "outcome_other_structural_abnormality")
  output$Supplemental_Fig_03 <- valve_rbind(lapply(secondary, function(outcome) {
    z <- valve_proportions(d, "Supplemental_Fig_03", outcome = outcome, include_all = TRUE, config = config)
    z$Denominator <- "All recorded TTE recipients in each mutually exclusive pathway"
    z
  }))
  output$Supplemental_Fig_04 <- valve_rbind(lapply(other, function(outcome) {
    z <- valve_proportions(d, "Supplemental_Fig_04", outcome = outcome, include_all = TRUE, config = config)
    z$Denominator <- "All recorded TTE recipients in each mutually exclusive pathway"
    z
  }))
  output$Ex_Data_Fig_07 <- valve_rbind(lapply(exam, function(outcome)
    valve_proportions(d, "Ex_Data_Fig_07", outcome = outcome, conditional = TRUE, include_all = TRUE, config = config)))
  output$Ex_Data_Fig_08 <- valve_rbind(lapply(management, function(outcome)
    valve_proportions(d, "Ex_Data_Fig_08", outcome = outcome, conditional = TRUE, include_all = TRUE, config = config)))
  for (i in 1:2) {
    figure <- sprintf("Supplemental_Fig_%02d", i)
    variables <- clinical[seq.int((i-1L)*4L+1L, i*4L)]
    if (i == 1L) variables <- c(condition, variables)
    output[[figure]] <- valve_proportions(d, figure, variables, config = config)
  }
  groups <- list(c(condition, clinical[1:4]), clinical[5:8])
  for (j in 1:2) for (g in 1:2) {
    figure <- sprintf("Supplemental_Fig_%02d", (j-1L)*2L+g+4L)
    outcome <- if (j == 1L) "outcome_vhd_composite_examinations" else "outcome_vhd_composite_managements"
    output[[figure]] <- valve_proportions(d, figure, groups[[g]], outcome, conditional = TRUE, config = config)
  }
  output
}

# VALVE-AI publication forest plots --------------------------------------------
VALVE_FOREST_OR_LIMITS <- c(0.25, 16)

.valve_forest_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  gsub("[\r\n]+", " ", x)
}

.valve_forest_number <- function(x, digits = 2L) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || is.na(x)) return("NE")
  if (is.infinite(x)) return("Inf")
  formatC(x, digits = digits, format = "f")
}

.valve_forest_p <- function(p, method = "") {
  p <- suppressWarnings(as.numeric(p))
  if (!length(p) || !is.finite(p)) return("NE")
  out <- if (p < 0.001) "<0.001" else sprintf("%.3f", p)
  if (grepl("Firth", .valve_forest_text(method), ignore.case = TRUE)) {
    out <- paste0(out, "\u2020")
  }
  out
}

.valve_forest_count <- function(e, n) {
  e <- suppressWarnings(as.numeric(e))
  n <- suppressWarnings(as.numeric(n))
  if (!length(n) || !is.finite(n) || n <= 0) return("NE")
  number <- function(z) format(round(z), big.mark = ",", scientific = FALSE,
                               trim = TRUE)
  if (!length(e) || !is.finite(e)) return(paste0("NA/", number(n)))
  paste0(number(e), "/", number(n), " (", sprintf("%.1f", 100 * e / n), "%)")
}

.valve_forest_layout <- function(df, figure_id) {
  required <- c("Outcome", "N_intervention", "Events_intervention", "N_control",
                "Events_control", "OR", "Lower", "Upper", "P")
  missing <- setdiff(required, names(df))
  if (length(missing)) stop("Forest source is missing: ", paste(missing, collapse = ", "))
  if (!nrow(df)) stop("No rows available for ", figure_id)
  for (key in c("Section", "Panel", "Population", "Subgroup", "Variable", "Level",
                "Outcome_variable", "Method", "Interaction_Method")) {
    if (!key %in% names(df)) df[[key]] <- ""
    df[[key]] <- .valve_forest_text(df[[key]])
  }
  if (!"Interaction_P" %in% names(df)) df$Interaction_P <- NA_real_
  if (!"Available" %in% names(df)) df$Available <- TRUE
  id <- sub("\\.(pdf|png)$", "", basename(figure_id), ignore.case = TRUE)
  id <- sub("^Fig_([2-5])$", "Fig_0\\1", id)
  supported <- c("Fig_02", "Fig_03", "Fig_04", "Fig_05",
                 "Ex_Data_Fig_02", "Ex_Data_Fig_03", "Ex_Data_Fig_04")
  if (!id %in% supported) stop("Not a forest figure identifier: ", id)
  interaction <- id %in% c("Fig_03", "Ex_Data_Fig_02", "Ex_Data_Fig_03", "Ex_Data_Fig_04")
  label_header <- if (id %in% c("Fig_02", "Fig_04", "Fig_05")) "Outcome" else
    if (id %in% c("Ex_Data_Fig_02", "Ex_Data_Fig_03")) "Subgroup" else "Outcome / subgroup"
  items <- list()
  add <- function(type, label, index = NA_integer_, shaded = FALSE) {
    items[[length(items) + 1L]] <<- list(type = type,
       label = .valve_forest_text(label), index = index, shaded = shaded)
  }
  indices <- seq_len(nrow(df))
  append_rows <- function(ii, label_column = "Outcome") {
    for (i in ii) add("row", df[[label_column]][i], i)
  }
  if (id == "Fig_02") {
    for (section in unique(df$Section)) {
      ii <- indices[df$Section == section]
      if (nzchar(section)) add("section", section, ii[1], TRUE)
      append_rows(ii)
    }
  } else if (id == "Fig_03") {
    # No a/b panels: utilization, primary and secondary outcomes are sections.
    buckets <- list("Utilization outcome" = df$Outcome_variable == "outcome_tte_after_randomization",
                    "Primary outcome" = df$Outcome_variable == "outcome_vhd",
                    "Secondary outcome" = !df$Outcome_variable %in%
                       c("outcome_tte_after_randomization", "outcome_vhd"))
    for (group in names(buckets)) {
      ii <- indices[buckets[[group]]]
      if (!length(ii)) next
      add("group", group, shaded = TRUE)
      for (outcome in unique(df$Outcome_variable[ii])) {
        jj <- ii[df$Outcome_variable[ii] == outcome]
        add("section", df$Outcome[jj[1]], jj[1])
        for (j in jj) {
          lab <- df$Population[j]
          if (!nzchar(lab)) lab <- paste("AI-ECG", tolower(df$Level[j]), "risk")
          add("row", lab, j)
        }
      }
    }
  } else if (id %in% c("Fig_04", "Fig_05")) {
    for (population in unique(df$Population)) {
      ii <- indices[df$Population == population]
      letter <- df$Panel[ii[1]]
      if (!nzchar(letter)) letter <- letters[match(population, unique(df$Population))]
      add("group", paste(letter, population), shaded = TRUE)
      for (section in unique(df$Section[ii])) {
        jj <- ii[df$Section[ii] == section]
        if (nzchar(section)) add("section", section, jj[1])
        append_rows(jj)
      }
    }
  } else if (id %in% c("Ex_Data_Fig_02", "Ex_Data_Fig_03")) {
    add("group", df$Population[1], shaded = TRUE)
    for (variable in unique(df$Variable)) {
      ii <- indices[df$Variable == variable]
      label <- df$Subgroup[ii[1]]
      if (!nzchar(label)) label <- df$Section[ii[1]]
      add("section", label, ii[1])
      append_rows(ii, "Level")
    }
  } else if (id == "Ex_Data_Fig_04") {
    for (variable in unique(df$Variable)) {
      ii <- indices[df$Variable == variable]
      letter <- letters[match(variable, unique(df$Variable))]
      label <- df$Subgroup[ii[1]]
      if (!nzchar(label)) label <- df$Section[ii[1]]
      add("group", paste(letter, label), shaded = TRUE)
      for (outcome in unique(df$Outcome_variable[ii])) {
        jj <- ii[df$Outcome_variable[ii] == outcome]
        add("section", df$Outcome[jj[1]], jj[1])
        append_rows(jj, "Level")
      }
    }
  }
  # 10.8-point data-row pitch, with restrained separation between main groups.
  y <- 0.305
  for (j in seq_along(items)) {
    if (items[[j]]$type == "group" && j > 1L) y <- y + 0.035
    items[[j]]$y <- y
    y <- y + if (items[[j]]$type == "group") 0.160 else 0.150
  }
  axis_y <- y + 0.015
  height <- axis_y + 0.300
  list(df = df, items = items, interaction = interaction,
       label_header = label_header, axis = axis_y, height = height, id = id)
}

valve_forest_dimensions <- function(df, figure_id) {
  # 183 mm journal double-column width. Plot content, rather than white space,
  # determines the height; a manuscript figure must not exceed 247 mm.
  layout <- .valve_forest_layout(df, figure_id)
  if (layout$height > 247 / 25.4) {
    stop("Forest has too many rows for a 247-mm page: ", figure_id,
         ". Split scientifically coherent sections; do not shrink all text.")
  }
  c(width = 183 / 25.4, height = layout$height)
}

valve_draw_forest <- function(df, figure_id, config = list()) {
  # Optional config$font_family (default sans); config$forest_cex multiplies
  # typography, while config$forest_label_cex changes left-column text only.
  # Scientific headings stay on ONE line even if longer than the label column.
  # The caller opens one PDF device for this figure.
  g <- .valve_forest_layout(df, figure_id)
  df <- g$df
  family <- if (!is.null(config$font_family)) config$font_family else "sans"
  mult <- if (!is.null(config$forest_cex)) config$forest_cex else 1
  label_mult <- if (!is.null(config$forest_label_cex)) config$forest_label_cex else 1
  palette <- c(estimate = "#3C5488", confidence = "#3C5488", section = "#F1F2F3",
               rule = "#6B7078", reference = "#A7ADB4")
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), ps = 10,
                family = family, xaxs = "i", yaxs = "i", xpd = NA,
                fg = "black", bg = "white")
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(g$height, 0))
  # Small outer gutters prevent endpoint markers and column labels being cut.
  x <- if (g$interaction) {
    c(label = 0.015, int = 0.349, ctrl = 0.456, plot1 = 0.529, plot2 = 0.671,
      orr = 0.752, p = 0.864, pi = 0.947)
  } else {
    c(label = 0.015, int = 0.388, ctrl = 0.505, plot1 = 0.578, plot2 = 0.754,
      orr = 0.843, p = 0.957, pi = 0.957)
  }
  limits <- VALVE_FOREST_OR_LIMITS
  log_limits <- log2(limits)
  tx <- function(z) x["plot1"] + (log2(z) - log_limits[1]) / diff(log_limits) *
    (x["plot2"] - x["plot1"])
  write <- function(xx, yy, label, pt = 6.2, bold = FALSE, left = FALSE,
                    maximum_width = NULL) {
    cex <- pt / 10 * mult
    if (!is.null(maximum_width) && is.character(label)) {
      span <- graphics::strwidth(label, cex = cex, font = if (bold) 2 else 1)
      if (is.finite(span) && span > maximum_width) cex <- cex * maximum_width / span
    }
    graphics::text(xx, yy, labels = label, adj = c(if (left) 0 else 0.5, 0.5),
                   cex = cex, font = if (bold) 2 else 1)
  }
  write(x["label"], 0.090, g$label_header, 6.0, TRUE, TRUE)
  write(x["int"], 0.090, "Intervention", 6.0, TRUE)
  write(x["ctrl"], 0.090, "Control", 6.0, TRUE)
  write(x["orr"], 0.090, "Odds ratio (95% CI)", 6.0, TRUE)
  write(x["p"], 0.090, expression(bolditalic(P)~bold(value)), 6.0)
  if (g$interaction) {
    write(x["pi"], 0.090, expression(bolditalic(P)~bold("for interaction")), 6.0)
  }
  graphics::segments(0.012, 0.190, 0.988, 0.190, col = palette["rule"], lwd = 0.7)
  graphics::segments(tx(1), 0.190, tx(1), g$axis,
                     col = palette["reference"], lwd = 0.7, lty = "22")
  for (item in g$items) {
    y <- item$y
    if (item$type != "row") {
      if (item$shaded) {
        graphics::rect(0.012, y - 0.071, 0.988, y + 0.071,
                       col = palette["section"], border = NA)
      }
      # Long variable headings are intentionally allowed across empty columns.
      write(x["label"] + if (item$type == "section") 0.006 else 0,
            y, item$label, 6.5 * label_mult, TRUE, TRUE)
      if (g$interaction && item$type == "section" && !is.na(item$index)) {
        i <- item$index
        write(x["pi"], y, .valve_forest_p(df$Interaction_P[i],
                                         df$Interaction_Method[i]), 6.0)
      }
      next
    }
    i <- item$index
    write(x["label"] + 0.010, y, item$label, 6.2 * label_mult, left = TRUE)
    available <- isTRUE(as.logical(df$Available[i]))
    if (!available) {
      write(x["int"], y, "Not available", 6.0)
      write(x["ctrl"], y, "Not available", 6.0)
      write(x["orr"], y, "Not available", 6.0)
      write(x["p"], y, "\u2014", 6.0)
      next
    }
    write(x["int"], y, .valve_forest_count(df$Events_intervention[i],
                                          df$N_intervention[i]), 6.0,
          maximum_width = if (g$interaction) 0.105 else 0.115)
    write(x["ctrl"], y, .valve_forest_count(df$Events_control[i], df$N_control[i]),
          6.0, maximum_width = if (g$interaction) 0.105 else 0.115)
    est <- suppressWarnings(as.numeric(df$OR[i]))
    lo <- suppressWarnings(as.numeric(df$Lower[i]))
    hi <- suppressWarnings(as.numeric(df$Upper[i]))
    both_zero <- isTRUE(as.numeric(df$Events_intervention[i]) == 0 &&
                       as.numeric(df$Events_control[i]) == 0)
    orci <- if (is.na(est) || both_zero) "NE" else paste0(
      .valve_forest_number(est), " (", .valve_forest_number(lo), "\u2013",
      .valve_forest_number(hi), ")")
    write(x["orr"], y, orci, 6.0,
          maximum_width = if (g$interaction) 0.155 else 0.165)
    write(x["p"], y, if (both_zero) "NE" else .valve_forest_p(df$P[i], df$Method[i]), 6.0)
    if (both_zero || is.na(est) || is.na(lo) || is.na(hi)) next
    if (hi < lo || lo < 0) next
    lower <- max(lo, limits[1])
    upper <- min(hi, limits[2])
    if (upper >= lower && is.finite(lower) && is.finite(upper)) {
      graphics::segments(tx(lower), y, tx(upper), y,
                         col = palette["confidence"], lwd = 1.0)
    }
    if (is.finite(est) && est >= limits[1] && est <= limits[2]) {
      graphics::points(tx(est), y, pch = 15, cex = 0.43,
                       col = palette["estimate"])
    }
    # Arrow tips are exactly on the log-axis boundary, centered at the CI y.
    # Both sides use the same physical geometry; the point estimate is not
    # clipped to the boundary, which would misrepresent an out-of-range OR.
    arrow_dx <- 0.006
    arrow_dy <- 0.018
    if (lo < limits[1]) {
      graphics::lines(c(x["plot1"] + arrow_dx, x["plot1"], x["plot1"] + arrow_dx),
                      c(y - arrow_dy, y, y + arrow_dy),
                      col = palette["confidence"], lwd = 1.0)
    }
    if (hi > limits[2]) {
      graphics::lines(c(x["plot2"] - arrow_dx, x["plot2"], x["plot2"] - arrow_dx),
                      c(y - arrow_dy, y, y + arrow_dy),
                      col = palette["confidence"], lwd = 1.0)
    }
  }
  graphics::segments(x["plot1"], g$axis, x["plot2"], g$axis,
                     col = palette["rule"], lwd = 0.6)
  ticks <- 2^seq(log_limits[1], log_limits[2])
  tick_labels <- vapply(ticks, function(z)
    sub("^0[.]", ".", format(z, trim = TRUE, scientific = FALSE)), character(1))
  for (j in seq_along(ticks)) {
    xx <- tx(ticks[j])
    graphics::segments(xx, g$axis, xx, g$axis + 0.030,
                       col = palette["rule"], lwd = 0.6)
    write(xx, g$axis + 0.085, tick_labels[j], 6.0)
  }
  write(mean(x[c("plot1", "plot2")]), g$axis + 0.205,
        "Intervention / control odds ratio", 6.0)
  invisible(g)
}

.vyg_palette <- function(config = list()) {
  default <- c("#3C5488", "#E64B35", "#00A087", "#8491B4", "#F39B7F")
  if (is.list(config) && length(config$pathway_colors) >= 5L)
    return(as.character(config$pathway_colors[1:5]))
  default
}
.vyg_paths <- function() c(
  "AI-ECG high risk | Intervention: trial-triggered",
  "AI-ECG high risk | Intervention: physician-triggered",
  "AI-ECG high risk | Control: physician-triggered",
  "AI-ECG low risk | Intervention: physician-triggered",
  "AI-ECG low risk | Control: physician-triggered")
.vyg_p <- function(p) {
  if (length(p) == 0L || !is.finite(p[1])) return("NE")
  if (p[1] < .001) "<0.001" else sprintf("%.3f", p[1])
}
.vyg_pd <- function(row) {
  if ("P_display" %in% names(row) && !is.na(row$P_display[1]) &&
      nzchar(as.character(row$P_display[1]))) return(as.character(row$P_display[1]))
  if ("Pathway_index" %in% names(row) && row$Pathway_index[1] == 1 && row$N[1] > 0)
    return("Ref.")
  .vyg_p(row$P)
}
.vyg_num <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
.vyg_text <- function(label, x, y, size = 6.5, face = "plain", col = "#16191C",
                      just = "left", ...) {
  grid::grid.text(label, x = x, y = y, just = just,
                  gp = grid::gpar(fontfamily = getOption("valve_ai.yield_font", "sans"), fontsize = size,
                                  fontface = face, col = col), ...)
}
.vyg_line <- function(x, y, col = "#626A73", lwd = .6, lty = 1) {
  grid::grid.lines(x = x, y = y, gp = grid::gpar(col = col, lwd = lwd, lty = lty))
}
.vyg_rect <- function(x, y, width, height, fill = "#EDF0F3", col = NA) {
  grid::grid.rect(x = x, y = y, width = width, height = height,
                  gp = grid::gpar(fill = fill, col = col))
}
.vyg_strip_suffix <- function(x) sub(" \\(high-risk participants only\\)$", "", as.character(x))
.vyg_p_label <- function(p) {
  shown <- .vyg_p(p)
  if (shown == "NE") return(expression(bolditalic(P)~"value = NE"))
  if (shown == "<0.001") return(expression(bolditalic(P)~"value <0.001"))
  bquote(bolditalic(P)~"value = "*.(shown))
}

valve_yield_dimensions <- function(df, id) {
  width <- 183 / 25.4
  if (id == "Ex_Data_Fig_06") {
    row_units <- nrow(df) + length(unique(df$Variable)) * (1.25 + .22)
    return(c(width = width, height = (3.1 * row_units + 10) / 25.4))
  }
  height_mm <- c(
    Ex_Data_Fig_05 = 167.0, Ex_Data_Fig_07 = 116.893, Ex_Data_Fig_08 = 133.764,
    Supplemental_Fig_01 = 212.0, Supplemental_Fig_02 = 188.122,
    Supplemental_Fig_03 = 100.022, Supplemental_Fig_04 = 133.764,
    Supplemental_Fig_05 = 212.0, Supplemental_Fig_06 = 188.122,
    Supplemental_Fig_07 = 212.0, Supplemental_Fig_08 = 188.122)
  if (!id %in% names(height_mm)) stop("No yield-figure dimensions for ", id)
  c(width = width, height = unname(height_mm[id]) / 25.4)
}

# Convert a physical vertical distance to the current viewport coordinates.
.vyg_mm <- function(mm) grid::convertHeight(grid::unit(mm, "mm"), "npc", valueOnly = TRUE)

.vyg_legend <- function(colors, y = .vyg_mm(3.5), n = 5L) {
  labels <- c("AI-ECG high risk\nIntervention: trial-triggered",
              "AI-ECG high risk\nIntervention: physician-triggered",
              "AI-ECG high risk\nControl: physician-triggered",
              "AI-ECG low risk\nIntervention: physician-triggered",
              "AI-ECG low risk\nControl: physician-triggered")
  for (i in seq_len(n)) {
    x <- (i - .5) / n
    grid::grid.points(x = x - .080, y = y + .002, pch = 16,
                      default.units = "npc", size = grid::unit(1.25, "mm"),
                      gp = grid::gpar(col = colors[i]))
    .vyg_text(labels[i], x - .066, y, size = 5.4,
              just = c("left", "center"))
  }
}

.vyg_items <- function(df, outcome = FALSE, letters = TRUE) {
  result <- list(); counter <- 0L
  add <- function(kind, label = "", row = NULL, weight = 1) {
    counter <<- counter + 1L
    result[[counter]] <<- list(kind = kind, label = label, row = row, weight = weight)
  }
  key <- if (outcome) "Outcome_variable" else "Variable"
  groups <- unique(as.character(df[[key]]))
  for (j in seq_along(groups)) {
    block <- df[as.character(df[[key]]) == groups[j], , drop = FALSE]
    label <- if (outcome) as.character(block$Outcome[1]) else as.character(block$Subgroup[1])
    if (letters && !outcome) label <- paste0(base::letters[j], "  ", label)
    add("heading", label, weight = 1.15)
    level_values <- if (outcome) "__ALL__" else unique(as.character(block$Level))
    for (level in level_values) {
      sub <- if (outcome) block else block[as.character(block$Level) == level, , drop = FALSE]
      if (!outcome) add("category", level, weight = .9)
      sub <- sub[order(sub$Pathway_index), , drop = FALSE]
      for (i in seq_len(nrow(sub))) add("row", row = sub[i, , drop = FALSE])
      add("space", weight = .2)
    }
    add("space", weight = .15)
  }
  result
}

.vyg_forest_body <- function(items, top_label = "Outcome / TTE pathway",
                              colors, axis_label = "Proportion (%)", size = 6.4,
                              title = NULL, subgroup_p = FALSE, compact = FALSE) {
  # Current viewport is an entire panel. Coordinate geometry is independent of
  # row values; intervals are plotted on the fixed probability scale [0,1].
  x0 <- .455; x1 <- .680; value_x <- .817; p_x <- .993
  header <- if (compact) 1 - .vyg_mm(if (is.null(title)) 1.9 else 5.5) else
    if (is.null(title)) .969 else .934
  title_y <- if (compact) 1 - .vyg_mm(1.5) else .986
  if (!is.null(title)) .vyg_text(title, .006, title_y, size = size + .3, face = "bold")
  if (nzchar(top_label)) .vyg_text(top_label, .006, header, size = size, face = "bold")
  .vyg_text(axis_label, (x0 + x1) / 2, header, size = size, face = "bold", just = "center")
  .vyg_text("n/N; % (95% CI)", value_x, header, size = size, face = "bold", just = "center")
  if (!subgroup_p)
    .vyg_text(expression(bolditalic(P)~"value"), p_x, header, size = size, just = "right")
  rule_y <- header - if (compact) .vyg_mm(2.3) else .028
  .vyg_line(c(.002, .998), rep(rule_y, 2), lwd = .75)
  bottom <- if (compact) .vyg_mm(3.4) else .048
  available <- rule_y - bottom - if (compact) .vyg_mm(1.4) else .022
  weights <- vapply(items, function(x) x$weight, numeric(1))
  scale <- available / sum(weights)
  if (compact && scale < .vyg_mm((size - .45) * 25.4 / 72 + .25))
    stop("Figure rows exceed the readable spacing of the manuscript layout.")
  centers <- rule_y - (if (compact) .vyg_mm(.7) else .012) -
    (cumsum(weights) - weights / 2) * scale
  for (tick in seq(0, 1, .25)) {
    x <- x0 + (x1 - x0) * tick
    .vyg_line(c(x, x), c(bottom, rule_y), col = if (tick == 0) "#8C959F" else "#D6DADF",
              lwd = .6, lty = 3)
    .vyg_line(c(x, x), c(bottom, bottom - if (compact) .vyg_mm(.7) else .008))
    .vyg_text(as.character(round(100 * tick)), x, bottom - if (compact) .vyg_mm(2.0) else .023,
              size = size - .4, just = "center")
  }
  .vyg_line(c(x0, x1), c(bottom, bottom), lwd = .75)
  for (j in seq_along(items)) {
    item <- items[[j]]; y <- centers[j]
    if (item$kind == "space") next
    if (item$kind == "heading") {
      .vyg_rect(.5, y, .998, scale * item$weight * .91)
      .vyg_text(.vyg_strip_suffix(item$label), .006, y, size = size, face = "bold")
      if (subgroup_p && !is.null(item$row))
        .vyg_text(.vyg_p_label(item$row$P), p_x, y, size = size, just = "right")
      next
    }
    if (item$kind == "category") {
      .vyg_text(item$label, .013, y, size = size - .2, face = "bold")
      next
    }
    row <- item$row
    index <- as.integer(row$Pathway_index[1]); if (is.na(index)) index <- 1L
    color <- colors[index]
    label <- if (subgroup_p) as.character(row$Level[1]) else as.character(row$Pathway[1])
    .vyg_text(label, .017, y, size = size - .45)
    if (is.finite(row$N[1]) && row$N[1] > 0 && is.finite(row$Proportion[1])) {
      lo <- pmax(0, pmin(1, row$Lower[1])); hi <- pmax(0, pmin(1, row$Upper[1]))
      if (is.finite(lo) && is.finite(hi)) {
        .vyg_line(x0 + (x1 - x0) * c(lo, hi), rep(y, 2), col = color, lwd = 1)
        for (lim in c(lo, hi)) .vyg_line(rep(x0 + (x1 - x0) * lim, 2),
                                        y + c(-.10, .10) * scale, col = color, lwd = .7)
      }
      grid::grid.points(x = x0 + (x1 - x0) * row$Proportion[1], y = y,
                        pch = 16, default.units = "npc", size = grid::unit(1.35, "mm"),
                        gp = grid::gpar(col = color))
      label <- sprintf("%s/%s; %.1f (%.1f-%.1f)", .vyg_num(row$Events[1]),
                       .vyg_num(row$N[1]), 100 * row$Proportion[1],
                       100 * row$Lower[1], 100 * row$Upper[1])
    } else label <- "0/0; NE"
    .vyg_text(label, value_x, y, size = size - .45, col = color, just = "center")
    if (!subgroup_p) .vyg_text(.vyg_pd(row), p_x, y, size = size - .45, just = "right")
  }
}

.vyg_standard <- function(df, id, config) {
  colors <- .vyg_palette(config)
  outcome <- id %in% c(sprintf("Ex_Data_Fig_%02d", 7:8), sprintf("Supplemental_Fig_%02d", 3:4))
  items <- .vyg_items(df, outcome = outcome, letters = !outcome)
  grid::pushViewport(grid::viewport(x = .5, y = grid::unit(7.2, "mm"),
    just = c("center", "bottom"), width = .976,
    height = grid::unit(1, "npc") - grid::unit(8.4, "mm")))
  .vyg_forest_body(items, top_label = if (outcome) "Outcome / TTE pathway" else "Subgroup / TTE pathway",
                   colors = colors, axis_label = if (id %in% sprintf("Supplemental_Fig_%02d", 1:4)) "Yield (%)" else "Proportion (%)",
                   size = if (nrow(df) > 40L) 6.05 else 6.5, compact = TRUE)
  grid::popViewport(); .vyg_legend(colors)
}

.vyg_trial_subgroups <- function(df, config) {
  if ("Variable_order" %in% names(df))
    df <- df[order(df$Variable_order, df$Level_order), , drop = FALSE]
  items <- list(); k <- 0L
  add <- function(x) {k <<- k + 1L; items[[k]] <<- x}
  for (variable in unique(as.character(df$Variable))) {
    block <- df[as.character(df$Variable) == variable, , drop = FALSE]
    add(list(kind = "heading", label = as.character(block$Subgroup[1]),
             row = block[1, , drop = FALSE], weight = 1.25))
    for (i in seq_len(nrow(block)))
      add(list(kind = "row", label = "", row = block[i, , drop = FALSE], weight = 1))
    add(list(kind = "space", label = "", row = NULL, weight = .22))
  }
  grid::pushViewport(grid::viewport(x = .5, y = grid::unit(1, "mm"),
    just = c("center", "bottom"), width = .976,
    height = grid::unit(1, "npc") - grid::unit(2, "mm")))
  .vyg_forest_body(items, top_label = "", colors = .vyg_palette(config),
                   axis_label = "Yield (%)", size = 7, subgroup_p = TRUE, compact = TRUE)
  grid::popViewport()
  # Intentionally no panel letters, Subgroup heading, pathway title or legend.
}

.vyg_detection_and_yield <- function(df, config) {
  tree <- df[as.character(df$Panel) == "a", , drop = FALSE]
  yields <- df[as.character(df$Panel) == "b", , drop = FALSE]
  stopifnot(nrow(tree) > 0L, nrow(yields) == 5L,
            all(yields$Variable == "All"))
  grid::pushViewport(grid::viewport(x = .5, y = grid::unit(1, "npc"),
    just = c("center", "top"), width = 1,
    height = grid::unit(119.992, "mm")))
  .vyg_tree(tree, config)
  .vyg_text("a", .006, 1 - .vyg_mm(1.5), size = 7, face = "bold")
  grid::popViewport()
  items <- list(list(kind = "heading", label = "All participants", row = NULL, weight = 1.12))
  yields <- yields[order(yields$Pathway_index), , drop = FALSE]
  for (i in seq_len(nrow(yields)))
    items[[length(items)+1L]] <- list(kind = "row", label = "", row = yields[i, , drop = FALSE], weight = 1)
  grid::pushViewport(grid::viewport(x = .5, y = grid::unit(7.2, "mm"),
    just = c("center", "bottom"), width = .976,
    height = grid::unit(38.0, "mm")))
  .vyg_forest_body(items, top_label = "TTE pathway", colors = .vyg_palette(config),
    axis_label = "Yield (%)", size = 6.3, title = "b", compact = TRUE)
  grid::popViewport()
  .vyg_legend(.vyg_palette(config))
}

.vyg_tree <- function(df, config) {
  pie_colors <- c("#4DBBD5", "#E64B35")
  branch_color <- "#D3D7DD"
  text_color <- .vyg_palette(config)[1]
  # Node centers and parent/child branch offsets are explicitly symmetric.
  # No geometry is scaled by event counts, avoiding asymmetric pie branches.
  branch <- function(x1, y1, x2, y2, n, parent_n) {
    t <- seq(0, 1, length.out = 81)
    x <- (1-t)^3*x1 + 3*(1-t)^2*t*(x1+.085) + 3*(1-t)*t^2*(x2-.085) + t^3*x2
    y <- (1-t)^3*y1 + 3*(1-t)^2*t*y1 + 3*(1-t)*t^2*y2 + t^3*y2
    .vyg_line(x, y, col = branch_color, lwd = .85)
    label <- if (parent_n > 0) sprintf("%.1f%%", 100*n/parent_n) else "NE"
    midx <- x[41]; midy <- y[41]
    grid::grid.rect(x = midx, y = midy, width = grid::stringWidth(label) + grid::unit(2, "mm"),
                    height = grid::unit(3.0, "mm"), gp = grid::gpar(fill = "white", col = NA))
    .vyg_text(label, midx, midy, size = 6.2, col = text_color, just = "center")
  }
  pie <- function(x, y, p, unverified = FALSE) {
    r <- grid::unit(2.1, "mm")
    grid::grid.circle(x = x, y = y, r = r,
                      gp = grid::gpar(fill = if (unverified) "#D3D7DD" else pie_colors[1],
                                      col = "white", lwd = .5))
    if (!unverified && is.finite(p) && p > 0) {
      theta <- seq(pi/2, pi/2 - 2*pi*min(p, 1), length.out = max(3, ceiling(100*p)))
      grid::grid.polygon(x = grid::unit(x, "npc") + r*c(0, cos(theta), 0),
                         y = grid::unit(y, "npc") + r*c(0, sin(theta), 0),
                         gp = grid::gpar(fill = pie_colors[2], col = "white", lwd = .5))
    }
  }
  node <- function(row, x, y, title) {
    source <- as.character(row$SOURCE[1])
    unknown <- grepl("unknown|missing", source, ignore.case = TRUE)
    no_tte <- grepl("no tte|without tte", source, ignore.case = TRUE)
    denominator <- if ("Observed" %in% names(row)) row$Observed[1] else row$Total[1]
    p <- if (denominator > 0) row$Events[1]/denominator else NA_real_
    pie(x, y, p, no_tte || unknown)
    if (no_tte || unknown) {
      label <- paste0(title, "\nn = ", .vyg_num(row$Total[1]), "; unverified")
    } else {
      ci <- if (is.finite(row$CI_lower[1]) && is.finite(row$CI_upper[1]))
        sprintf("95%% CI %.1f%%-%.1f%%", 100*row$CI_lower[1], 100*row$CI_upper[1]) else "95% CI NE"
      label <- paste0(title, "\n", .vyg_num(row$Events[1]), "/", .vyg_num(denominator),
                      "; ", if (is.finite(p)) sprintf("%.1f%%", 100*p) else "NE", "\n", ci)
    }
    .vyg_text(label, x, y - .vyg_mm(3.0), size = 6.0, just = c("center", "top"))
  }
  for (arm in c("Intervention", "Control")) {
    offset <- if (arm == "Intervention") 0 else .5
    grid::pushViewport(grid::viewport(x = offset + .25, y = .53, width = .49, height = .90))
    .vyg_text(arm, .5, .985, size = 8.5, face = "bold", just = "center")
    block <- df[as.character(df$randomization_group) == arm, , drop = FALSE]
    all <- block[block$Layer == 1, , drop = FALSE]
    if (!nrow(all)) {grid::popViewport(); next}
    root_x <- .11; middle_x <- .39; end_x <- .81; root_y <- .50
    risk_y <- c(High = .73, Low = .27)
    risks <- block[block$Layer == 2, , drop = FALSE]
    positions <- list()
    # Draw all branches first so labels and pies stay above the connectors.
    for (risk in names(risk_y)) {
      rr <- risks[tolower(as.character(risks$AI_RISK)) == tolower(risk), , drop = FALSE]
      if (!nrow(rr)) next
      yy <- risk_y[[risk]]
      branch(root_x, root_y, middle_x, yy, rr$Total[1], all$Total[1])
      children <- block[block$Layer == 3 & tolower(as.character(block$AI_RISK)) == tolower(risk), , drop = FALSE]
      source <- as.character(children$SOURCE)
      ord <- order(match(tolower(source), c("trial-triggered", "physician-triggered", "no tte", "unknown")), na.last = TRUE)
      children <- children[ord, , drop = FALSE]
      n <- nrow(children)
      offsets <- if (n <= 1) 0 else seq(.112, -.112, length.out = n)
      if (n == 3) offsets <- c(.145, 0, -.145)
      if (n >= 4) offsets <- seq(.16, -.16, length.out = n)
      for (j in seq_len(n)) {
        cy <- yy + offsets[j]
        branch(middle_x, yy, end_x, cy, children$Total[j], rr$Total[1])
        positions[[length(positions)+1L]] <- list(row = children[j, , drop = FALSE], y = cy)
      }
    }
    node(all[1, , drop = FALSE], root_x, root_y, "All participants")
    for (risk in names(risk_y)) {
      rr <- risks[tolower(as.character(risks$AI_RISK)) == tolower(risk), , drop = FALSE]
      if (nrow(rr)) node(rr[1, , drop = FALSE], middle_x, risk_y[[risk]], paste("AI-ECG", tolower(risk), "risk"))
    }
    for (pos in positions) {
      src <- as.character(pos$row$SOURCE[1])
      title <- if (tolower(src) == "trial-triggered") "Trial-triggered TTE" else
        if (tolower(src) == "physician-triggered") "Physician-triggered TTE" else
          if (tolower(src) == "no tte") "Without TTE" else "Unknown TTE status"
      node(pos$row, end_x, pos$y, title)
    }
    grid::popViewport()
  }
  labels <- c("No primary detection recorded", "Primary detection", "No TTE: disease status unverified")
  cols <- c(pie_colors, "#D3D7DD")
  centers <- c(.25, .52, .80)
  for (i in 1:3) {
    .vyg_rect(centers[i]-.13, .033, .022, .012, fill = cols[i])
    .vyg_text(labels[i], centers[i]-.112, .033, size = 6.1)
  }
}

# Main drawing entry point. The caller opens a PDF device using
# valve_yield_dimensions() and then calls this function exactly once per page.
valve_draw_yield <- function(df, figure_id, config = list()) {
  if (!is.data.frame(df) || !nrow(df)) stop("No source rows for ", figure_id)
  family <- if (is.list(config) && !is.null(config$font_family)) config$font_family else "sans"
  old_options <- options(valve_ai.yield_font = family)
  on.exit(options(old_options), add = TRUE)
  grid::grid.newpage()
  if (figure_id == "Ex_Data_Fig_05") return(invisible(.vyg_detection_and_yield(df, config)))
  if (figure_id == "Ex_Data_Fig_06") return(invisible(.vyg_trial_subgroups(df, config)))
  invisible(.vyg_standard(df, figure_id, config))
}

# Descriptive ECG timing and frequency panels; no hypothesis tests.
valve_ecg_source <- function(d) {
  # Integer-day bins: 0-15, then 16-30, 31-45, ..., 361-375.
  # The first bin exactly matches the 0-15-day inset; eligible days end at 365.
  out <- list(); breaks <- c(0, seq(16,376,15))
  for (panel in c("a","b","c","d")) for (arm in c("Intervention","Control")) {
    mask <- d$randomization_group==arm
    if (panel=="b") mask <- mask & d$AI_HIGH
    values <- switch(panel,a=d$days_since_last_ecg[mask],b=d$days_since_last_high_risk_ecg[mask],
                     c=unlist(attr(d,"ecg_lists")[mask],use.names=FALSE),d=d$n_ECG[mask])
    values <- values[!is.na(values)]
    qu <- if(length(values)) as.numeric(quantile(values,c(.5,.25,.75),type=7)) else rep(NA_real_,3)
    if (panel!="d") {
      if(any(values<0 | values>365)) stop("ECG day intervals outside 0-365 in panel ",panel)
      counts <- hist(values,breaks=breaks,plot=FALSE,right=FALSE,include.lowest=TRUE)$counts
      tab <- data.frame(Bin_lower=head(breaks,-1),Bin_upper=tail(breaks,-1),Count=counts,Category=NA_character_,Category_order=NA_integer_)
      rule <- "Integer days [lower, upper): 0-15, 16-30, 31-45, ...; eligible days 0-365"
    } else {
      counts <- tabulate(pmin(values,10),nbins=10)
      tab <- data.frame(Bin_lower=NA_real_,Bin_upper=NA_real_,Count=counts,Category=c(as.character(1:9),">=10"),Category_order=1:10)
      rule <- "Exact count; final category >=10"
    }
    tab$Figure <- "Ex_Data_Fig_01";tab$Panel <- panel;tab$Group <- arm
    tab$Variable <- switch(panel,a="DAYS_RECENT_ECG",b="DAYS_RECENT_HIGH_RISK_ECG",c="Days from all available ECGs",d="Number of available ECGs per participant")
    tab$Unit <- if(panel=="c")"ECGs" else "Participants";tab$N <- length(values)
    tab$Median <- qu[1];tab$Q1<-qu[2];tab$Q3<-qu[3];tab$Proportion<-if(length(values))counts/length(values) else NA_real_
    tab$Bin_rule<-rule;tab$Input_column<-switch(panel,a="days_since_last_ecg",b="days_since_last_high_risk_ecg",c="ecg_days_before_index",d="count of entries in ecg_days_before_index")
    tab$View <- "Overview"
    out[[length(out)+1L]] <- tab
    if (panel != "d") {
      # Zoom uses integer-day bins 0-1, 2-3, ..., 14-15. Day 16 is
      # excluded; percentages retain the full parent-panel denominator.
      if (any(values != floor(values))) stop("Non-integer ECG days in panel ", panel)
      zoom <- tab[seq_len(8L), , drop = FALSE]
      zoom$Bin_lower <- seq(0, 14, 2)
      zoom$Bin_upper <- zoom$Bin_lower + 2
      zoom$Count <- vapply(zoom$Bin_lower, function(lo)
        sum(values >= lo & values < lo + 2), integer(1))
      zoom$Category <- paste0(zoom$Bin_lower, "-", zoom$Bin_upper - 1)
      zoom$Category_order <- seq_len(8L)
      zoom$Proportion <- if (length(values)) zoom$Count / length(values) else NA_real_
      zoom$Bin_rule <- "Integer days [lower, upper); 0-15 included, day 16 excluded"
      zoom$View <- "Inset_0_15_days"
      out[[length(out)+1L]] <- zoom
    }
  }
  do.call(rbind,out)
}

valve_draw_ecg <- function(df, config) {
  # Counts and each legend denominator come only from valve_ecg_source().
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(1, 1), mar = rep(0, 4), oma = rep(0, 4),
    family = config$font_family, ps = 10, cex = 1, xaxs = "i", yaxs = "i",
    mgp = c(3, .65, 0), tcl = -.22, las = 1, fg = "#363435", bg = "white")
  colors <- c(Intervention = "#485E8F", Control = "#E85543")
  regions <- list(a = c(.0776, .4741, .5878, .9403),
                  b = c(.5884, .9850, .5878, .9403),
                  c = c(.0776, .4741, .1170, .4379),
                  d = c(.5884, .9850, .1170, .4379))
  xlabels <- c(a = "Days from most recent ECG",
               b = "Days from most recent high-risk ECG",
               c = "Days from all available ECGs",
               d = "Number of ECGs per participant")
  steps <- c(a = 200, b = 20, c = 250, d = 500)
  page_text <- function(label, x, y, size, font = 1, srt = 0, adj = .5) {
    usr <- graphics::par("usr"); area <- graphics::par("plt")
    xx <- usr[1] + (x - area[1]) / diff(area[1:2]) * diff(usr[1:2])
    yy <- usr[3] + (y - area[3]) / diff(area[3:4]) * diff(usr[3:4])
    graphics::text(xx, yy, label, cex = size / 10, font = font,
      srt = srt, adj = adj, xpd = NA, col = "#363435")
  }
  for (panel in c("a", "b", "c", "d")) {
    q <- df[df$Panel == panel & df$View == "Overview", , drop = FALSE]
    ymax <- max(q$Count, 1) * 1.23
    yticks <- seq(0, floor(ymax / steps[panel]) * steps[panel], by = steps[panel])
    graphics::par(plt = regions[[panel]], new = panel != "a")
    graphics::plot.new()
    graphics::plot.window(xlim = if (panel == "d") c(.4, 10.6) else c(-3, 378),
      ylim = c(0, ymax), xaxs = "i", yaxs = "i")
    if (panel != "d") {
      # Shade the first bin and its expansion behind all bars and percentages.
      # Connect vertically ordered points on the bin's right edge to the
      # corresponding inset endpoints, avoiding crossings and bar labels.
      first <- q[q$Bin_lower == 0, , drop = FALSE]
      highlight_top <- max(first$Count) + ymax * .13
      highlight_anchor <- max(first$Count) * .90
      graphics::rect(0, 0, 16, highlight_top, col = "#F0F1F3", border = NA)
      guide_x <- -3 + (.28 + .07) * 381
      guide_y <- ymax * c(1 - .02 - .12, .26 + .085)
      graphics::polygon(c(16,16,guide_x,guide_x),
        c(highlight_top,highlight_anchor,rev(guide_y)), col = "#F0F1F3", border = NA)
      graphics::segments(c(16,16), c(highlight_top,highlight_anchor), rep(guide_x,2),
        guide_y, col = "#A9AFB8", lwd = .55)
    }
    for (arm in names(colors)) {
      z <- q[q$Group == arm, , drop = FALSE]
      if (panel == "d") {
        left <- z$Category_order + if (arm == "Intervention") -.34 else 0
        right <- left + .34
      } else {
        midpoint <- (z$Bin_lower + z$Bin_upper) / 2
        left <- if (arm == "Intervention") z$Bin_lower + 1.1 else midpoint
        right <- if (arm == "Intervention") midpoint else z$Bin_upper - 1.1
      }
      graphics::rect(left, 0, right, z$Count, col = colors[arm], border = NA)
      # Direct group-specific percentages; vertical labels fit the paired
      # overview bars without merging neighboring bins. Preserve all zero bins.
      label_x <- (left + right) / 2
      graphics::text(label_x, z$Count + ymax * .018,
        labels = sprintf("%.1f%%", 100 * z$Proportion),
        srt = 90, adj = c(0, .5), cex = if (panel == "d") .486 else .405,
        col = colors[arm], xpd = NA)
    }
    xticks <- if (panel == "d") 1:10 else seq(0, 360, 60)
    xvalues <- if (panel == "d") c(as.character(1:9), "\u226510") else xticks
    graphics::axis(1, at = xticks, labels = FALSE, col = "#363435", lwd = .8)
    graphics::axis(2, at = yticks, labels = FALSE, col = "#363435", lwd = .8)
    usr <- graphics::par("usr"); area <- regions[[panel]]
    tick_x <- area[1] + (xticks - usr[1]) / diff(usr[1:2]) * diff(area[1:2])
    tick_y <- area[3] + (yticks - usr[3]) / diff(usr[3:4]) * diff(area[3:4])
    page_text(xvalues, tick_x, area[3] - 2.7 / 124.037, size = 6.2)
    page_text(yticks, area[1] - 2 / 183, tick_y, size = 6.2, adj = c(1, .5))
    graphics::box(bty = "l", col = "#363435", lwd = .8)
    counts <- vapply(names(colors), function(arm) {
      n <- unique(q$N[q$Group == arm])
      if (length(n) != 1L) stop("Inconsistent ECG legend denominator in panel ", panel)
      n
    }, numeric(1))
    legend_labels <- paste0(names(colors), " (n = ",
      format(counts, big.mark = ",", scientific = FALSE, trim = TRUE),
      if (panel == "c") " ECGs)" else ")")
    left_column <- panel %in% c("a", "c")
    top_row <- panel %in% c("a", "b")
    # Put both arm labels side by side below each panel's x-axis title.
    legend_page_y <- if (top_row) .495 else .020
    legend_user_y <- usr[3] + (legend_page_y-area[3])/diff(area[3:4])*diff(usr[3:4])
    graphics::legend(x = mean(usr[1:2]), y = legend_user_y,
      xjust = .5, yjust = .5, horiz = TRUE, xpd = NA,
      legend = legend_labels, fill = colors, border = NA, bty = "n",
      cex = .58, text.col = "#363435", x.intersp = .65, seg.len = 1.1)
    page_text(panel, if (left_column) .0348 else .5456,
      if (top_row) .9667 else .4645, size = 6.8, font = 2, adj = 0)
    page_text(xlabels[panel], mean(regions[[panel]][1:2]),
      if (top_row) .535 else .060, size = 7)
    page_text(if (panel == "c") "ECGs (n)" else "Participants (n)",
      if (left_column) .0285 else .5393, mean(regions[[panel]][3:4]), size = 7, srt = 90)
    if (panel != "d") {
      # The inset occupies the empty upper-right ~70% of the parent plot.
      # Keep its lower edge above all underlying bars and percentage labels.
      z <- df[df$Panel == panel & df$View == "Inset_0_15_days", , drop = FALSE]
      area <- regions[[panel]]
      xspan <- diff(area[1:2]); yspan <- diff(area[3:4])
      outer <- c(area[1] + .28 * xspan, area[2] - .01 * xspan,
                 area[3] + .26 * yspan, area[4] - .02 * yspan)
      usr <- graphics::par("usr")
      ux <- function(x) usr[1] + (x-area[1])/xspan*diff(usr[1:2])
      uy <- function(y) usr[3] + (y-area[3])/yspan*diff(usr[3:4])
      page_text("0-15 days", mean(outer[1:2]),
        outer[4] - .045*yspan, size = 5.7)
      inner <- c(outer[1] + .07*xspan, outer[2] - .02*xspan,
                 outer[3] + .085*yspan, outer[4] - .12*yspan)
      graphics::par(plt = inner, new = TRUE, tcl = -.13)
      graphics::plot.new()
      zymax <- max(z$Count, 1) * 1.29
      graphics::plot.window(xlim = c(-.2,16.2), ylim = c(0,zymax), xaxs = "i", yaxs = "i")
      zyticks <- pretty(c(0,zymax), n = 3)
      zyticks <- zyticks[zyticks >= 0 & zyticks <= zymax]
      for (arm in names(colors)) {
        zz <- z[z$Group == arm, , drop = FALSE]
        left <- zz$Bin_lower + if (arm == "Intervention") .18 else 1
        right <- left + .82
        graphics::rect(left, 0, right, zz$Count, col = colors[arm], border = NA)
        graphics::text((left+right)/2, zz$Count + zymax*.025,
          sprintf("%.1f%%",100*zz$Proportion), srt = 90, adj = c(0,.5),
          cex = .405, col = colors[arm], xpd = NA)
      }
      graphics::axis(1, at = c(0,5,10,15), labels = FALSE, lwd = .6, col = "#363435")
      graphics::axis(2, at = zyticks, labels = FALSE, lwd = .6, col = "#363435")
      graphics::box(bty = "l", lwd = .6, col = "#363435")
      page_text(c(0,5,10,15),
        inner[1] + (c(0,5,10,15)+.2)/16.4*diff(inner[1:2]),
        inner[3] - .035*yspan, size = 4.8)
      page_text(zyticks, inner[1] - .012*xspan,
        inner[3] + zyticks/zymax*diff(inner[3:4]), size = 4.8, adj = c(1,.5))
      graphics::par(tcl = -.22)
    }
  }
  invisible(NULL)
}

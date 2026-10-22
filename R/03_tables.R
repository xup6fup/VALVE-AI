
valve_table_p <- function(p, chi_square_fallback=FALSE) {
  if (!length(p) || is.na(p) || !is.finite(p)) return("NE")
  s <- if (p <= 0) paste0("<", format(.Machine$double.xmin, scientific=TRUE)) else if (p < .001) {
    z <- strsplit(sprintf("%.2e", p), "e", fixed=TRUE)[[1]]
    paste0(z[1], "×10^", as.integer(z[2]))
  } else sprintf("%.3f", p)
  paste0(s, if (chi_square_fallback) "§" else "")
}

valve_table_fisher <- function(tab, config = list()) {
  workspace <- if (is.null(config$table_fisher_workspace)) 200000L else config$table_fisher_workspace
  tab <- as.matrix(tab)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  ans <- list(P_value = NA_real_, Chi_square_fallback = FALSE)
  if (nrow(tab) < 2L) return(ans)
  if (ncol(tab) < 2L) {
    # Retain a zero-count complement for constant categorical variables.
    tab <- cbind(tab, `Zero-count complement` = 0)
  }
  exact <- tryCatch(suppressWarnings(
    stats::fisher.test(tab, alternative = "two.sided", workspace = workspace,
                       simulate.p.value = FALSE)), error = function(e) NULL)
  if (!is.null(exact)) {
    ans$P_value <- exact$p.value
    return(ans)
  }
  # Use Pearson's test only after an exact Fisher calculation fails.
  ans$Chi_square_fallback <- TRUE
  pearson <- tryCatch(suppressWarnings(
    stats::chisq.test(tab, correct = TRUE, simulate.p.value = FALSE)),
    error = function(e) NULL)
  if (!is.null(pearson)) ans$P_value <- pearson$p.value
  ans
}

valve_table_baseline <- function(d) {
  # Use the supplied clinical flags: 1 denotes a recorded condition,
  # 0 denotes no recorded condition, and missing values remain missing.
  flag <- function(col) ifelse(is.na(d[[col]]), NA_integer_,
                              as.integer(d[[col]] %in% c(1, "1")))
  map <- c(
    aneurysm = "history_aortic_aneurysm", # Recorded aortic aneurysm.
    marfan = "history_marfan_syndrome",   # Recorded Marfan syndrome.
    ef40 = "history_lvef_le40",     # Prior documented LVEF <=40%.
    af_recent = "atrial_fibrillation_within_3y",# AF with recent-onset (<=3-year) classification.
    hf_recent = "heart_failure_within_3y",# HF with recent-onset (<=3-year) classification.
    ckd_recent = "ckd_within_3y",# CKD with recent-onset (<=3-year) classification.
    cad_recent = "coronary_disease_within_3y",# CAD with recent-onset (<=3-year) classification.
    mi_recent = "myocardial_infarction_within_3y",# MI with recent-onset (<=3-year) classification.
    nt_high = "ntprobnp_ge900_within_3y",  # NT-proBNP >=900 pg/mL within 3 years.
    dyspnea = "dyspnea_within_3m",  # Dyspnea documented in preceding 3 months.
    palpitations = "palpitations_within_3m", # Palpitations in preceding 3 months.
    angina = "angina_chest_pain_within_3m",   # Angina/chest pain in preceding 3 months.
    syncope = "syncope_within_3m",  # Syncope in preceding 3 months.
    orthopnea = "orthopnea_within_3m",# Orthopnea in preceding 3 months.
    edema = "edema_within_3m",    # Edema in preceding 3 months.
    dm75 = "diabetes_age_ge75",        # Supplied age >=75 plus DM indicator.
    htn75 = "hypertension_age_ge75",       # Supplied age >=75 plus HTN indicator.
    ckd75 = "ckd_over_3y_age_ge75",       # Supplied age >=75 plus CKD indicator.
    thyroid = "history_thyroid_disorder",     # Recorded thyroid disorder.
    af_long = "atrial_fibrillation_over_3y",     # Long-standing AF (>3 years).
    hf_long = "heart_failure_over_3y",     # Long-standing HF (>3 years).
    cad_long = "coronary_disease_over_3y",    # Long-standing CAD (>3 years).
    mi_long = "myocardial_infarction_over_3y",     # Long-standing MI (>3 years).
    pacemaker = "history_pacemaker",   # Recorded permanent pacemaker.
    nt_moderate = "history_elevated_ntprobnp", # NT-proBNP >=300 in sinus rhythm / >=600 in AF.
    pasp = "history_elevated_pasp",        # Supplied elevated PASP history indicator.
    septum = "history_lv_hypertrophy",      # Supplied marked septal hypertrophy indicator.
    effusion = "history_large_pericardial_effusion",    # Supplied large pericardial effusion indicator.
    dm = "history_diabetes",         # Supplied DM history, independent of age.
    htn = "history_hypertension",       # Supplied HTN history, independent of age.
    ckd_long = "ckd_over_3y",  # CKD history outside recent-onset category.
    other = "history_other_structural_abnormality",           # Other structural/pericardial cardiac history.
    psych = "psychiatry_visit_within_180d",           # Psychiatric consultation in prior 6 months.
    ace = "baseline_acei_arb",             # Qualifying ACE inhibitor or ARB prescription.
    arni = "baseline_arni",            # Qualifying ARNI prescription.
    bb = "baseline_beta_blocker",              # Qualifying beta blocker prescription.
    mra = "baseline_mra",             # Qualifying mineralocorticoid antagonist.
    sglt2 = "baseline_sglt2_inhibitor",           # Qualifying SGLT2 inhibitor prescription.
    loop = "baseline_loop_diuretic",            # Qualifying loop diuretic prescription.
    ccb = "baseline_calcium_channel_blocker",             # Qualifying calcium channel blocker.
    thiazide = "baseline_thiazide",        # Qualifying thiazide prescription.
    alpha = "baseline_alpha_blocker",           # Qualifying alpha blocker prescription.
    otherhtn = "baseline_other_antihypertensive",        # Other qualifying antihypertensive drug.
    cv = "regular_cardiology_followup",              # Regular cardiology follow-up indicator.
    cvs = "regular_cardiovascular_surgery_followup"              # Regular cardiovascular-surgery follow-up.
  )
  need <- c(unname(map), "sex", "age_years", "age_group",
            "ai_ecg_risk", "days_since_last_ecg", "days_since_last_high_risk_ecg", "history_lvef_below50",
            "clinical_condition", c("vhd_valve_count_category", "vhd_aortic_stenosis_severity", "vhd_aortic_regurgitation_severity", "vhd_mitral_regurgitation_severity", "vhd_tricuspid_regurgitation_severity", "vhd_lvef_percent", "vhd_lvef_category", "vhd_lvedd_mm", "vhd_lvedd_category", "vhd_pasp_mmhg", "vhd_pasp_category", "vhd_postindex_dyspnea", "vhd_postindex_palpitations", "vhd_postindex_angina_chest_pain", "vhd_postindex_syncope", "vhd_postindex_orthopnea", "vhd_postindex_edema"))
  missing <- setdiff(need, names(d))
  if (length(missing)) stop("Missing table variables: ", paste(missing, collapse = ", "))
  x <- as.data.frame(lapply(map, flag), stringsAsFactors = FALSE)
  # Read the supplied demographic variables and age groups (60-74 / 75-85).
  x$male <- as.integer(tolower(as.character(d$sex)) == "male")
  x$age <- as.numeric(d$age_years)
  x$age_group <- as.character(d$age_group)
  x$age60 <- as.integer(x$age_group == "60-74")
  x$age75 <- as.integer(x$age_group == "75-85")
  x$high <- as.integer(d$ai_ecg_risk == "High")
  x$low <- as.integer(d$ai_ecg_risk == "Low")
  # days_since_last_ecg = days from most recent ECG to index time; days_since_last_high_risk_ecg = days from most
  # recent high-risk ECG, with NA in participants with no high-risk recording.
  x$days_ecg <- as.numeric(d$days_since_last_ecg)
  x$days_high <- as.numeric(d$days_since_last_high_risk_ecg)
  # history_lvef_below50 includes any supplied LVEF <50% history. Subtract the <=40%
  # category to obtain the mutually exclusive 41-49% category used in Table 1.
  x$ef41 <- as.integer(flag("history_lvef_below50") == 1 & x$ef40 == 0)
  x$ckd <- as.integer(x$ckd_recent == 1 | x$ckd_long == 1)
  x$condition <- as.character(d$clinical_condition)
  conditions <- c("High-priority echocardiography condition",
                  "Potential echocardiography condition",
                  "No prespecified echocardiography condition")
  for (i in seq_along(conditions)) x[[paste0("condition", i)]] <-
    as.integer(x$condition == conditions[i])
  # Recent, long-standing and unrecorded history form ONE global three-level
  # categorical test. The displayed rows show the first two categories only.
  for (nm in c("af", "hf", "ckd", "cad", "mi")) {
    recent <- x[[paste0(nm, "_recent")]] == 1
    longstanding <- x[[paste0(nm, "_long")]] == 1
    x[[paste0(nm, "_group")]] <- ifelse(is.na(recent) | is.na(longstanding),
      NA_character_, ifelse(recent, "Recent", ifelse(longstanding,
        "Long-standing", "No recorded history")))
  }
  # 'No documented reduction' must not be interpreted as measured LVEF >=50%.
  x$ef_group <- ifelse(x$ef40 == 1, "<=40%", ifelse(x$ef41 == 1, "41-49%",
                                                 "No documented reduction"))
  # vhd_valve_count_category: supplied count category of moderate/severe valve lesions.
  # Severity, echocardiographic measurements and categories are used as supplied.
  # Post-index symptom flags require documentation after VHD diagnosis and
  # within 90 days after index time; these final flags are used directly.
  for (nm in c("vhd_valve_count_category", "vhd_aortic_stenosis_severity", "vhd_aortic_regurgitation_severity", "vhd_mitral_regurgitation_severity", "vhd_tricuspid_regurgitation_severity", "vhd_lvef_percent", "vhd_lvef_category", "vhd_lvedd_mm", "vhd_lvedd_category", "vhd_pasp_mmhg", "vhd_pasp_category", "vhd_postindex_dyspnea", "vhd_postindex_palpitations", "vhd_postindex_angina_chest_pain", "vhd_postindex_syncope", "vhd_postindex_orthopnea", "vhd_postindex_edema")) x[[nm]] <- d[[nm]]
  x
}

valve_table_summary <- function(v, kind = "binary", category = 1) {
  v <- v[!is.na(v)]
  if (!length(v)) return("NE")
  if (kind == "mean") {
    v <- as.numeric(v)
    sd <- stats::sd(v)
    return(sprintf("%.1f ± %s", mean(v),
                    if (is.na(sd)) "NE" else sprintf("%.1f", sd)))
  }
  if (kind == "median") {
    q <- stats::quantile(as.numeric(v), c(.25, .5, .75), type = 7, names = FALSE)
    # Integer ECG intervals: arithmetic half-up rounding for nonnegative days.
    r <- floor(q + .5)
    return(sprintf("%.0f (%.0f-%.0f)", r[2], r[1], r[3]))
  }
  events <- sum(as.character(v) == as.character(category))
  paste0(format(events, big.mark = ",", scientific = FALSE, trim = TRUE),
         " (", sprintf("%.1f", 100 * events / length(v)), ")")
}

valve_table_test <- function(x, groups, key, kind, config = list()) {
  vv <- lapply(groups, function(g) {v <- x[[key]][g$mask]; v[!is.na(v)]})
  vv <- vv[lengths(vv) > 0L]
  rec <- list(P_value = NA_real_, Chi_square_fallback = FALSE)
  if (length(vv) < 2L) return(rec)
  if (kind %in% c("student", "anova")) {
    z <- tryCatch({
      if (kind == "student") {
        if (length(vv) != 2L) stop("Student t test needs two observed groups")
        stats::t.test(as.numeric(vv[[1]]), as.numeric(vv[[2]]),
                      var.equal = TRUE, alternative = "two.sided")$p.value
      } else {
        y <- as.numeric(unlist(vv, use.names = FALSE))
        group <- factor(rep(seq_along(vv), lengths(vv)))
        stats::anova(stats::lm(y ~ group))[["Pr(>F)"]][1]
      }
    }, error = function(e) NA_real_)
    rec$P_value <- if (is.finite(z)) as.numeric(z) else NA_real_
    return(rec)
  }
  cats <- sort(unique(as.character(unlist(vv, use.names = FALSE))))
  # Keep an explicit group-by-category matrix when a variable is constant.
  # vapply(..., numeric(1)) would simplify to a vector and transpose the wrong
  # dimension, incorrectly labelling all-zero baseline variables as NE.
  tab <- do.call(rbind, lapply(vv, function(v)
    as.numeric(table(factor(as.character(v), levels = cats)))))
  valve_table_fisher(tab, config)
}

valve_tables <- function(d, config = list()) {
  # AI_HIGH and PATHWAY_INDEX identify analysis groups from the prepared data.
  # PATHWAY_INDEX: 1 high/I/trial, 2 high/I/physician, 3 high/C/physician,
  # 4 low/I/physician, 5 low/C/physician; only recorded TTE recipients have a path.
  need <- c("randomization_group", "AI_HIGH", "PATHWAY_INDEX", "outcome_vhd",
            "outcome_tte", "tte_trigger")
  if (length(setdiff(need, names(d)))) stop("Run the input preparation module before valve_tables().")
  x <- valve_table_baseline(d)
  spec <- valve_table_schema()
  conditions <- c("High-priority echocardiography condition",
                  "Potential echocardiography condition",
                  "No prespecified echocardiography condition")
  paths <- c("AI-ECG high risk | Intervention: trial-triggered",
             "AI-ECG high risk | Intervention: physician-triggered",
             "AI-ECG high risk | Control: physician-triggered",
             "AI-ECG low risk | Intervention: physician-triggered",
             "AI-ECG low risk | Control: physician-triggered")
  group <- function(label, mask, stratum = "") list(label = label, mask = mask %in% TRUE,
                                                   stratum = stratum)
  gmain <- list(); gs1 <- list()
  for (risk in c("All participants", "High risk identified by AI-ECG", "Low risk identified by AI-ECG")) {
    keep <- if (risk == "All participants") rep(TRUE, nrow(d)) else if
      (risk == "High risk identified by AI-ECG") d$AI_HIGH == 1 else d$AI_HIGH == 0
    for (arm in c("Intervention", "Control")) gmain[[length(gmain) + 1L]] <-
      group(paste(risk, arm, sep = " | "), keep & d$randomization_group == arm, risk)
  }
  for (cc in conditions) for (arm in c("Intervention", "Control"))
    gs1[[length(gs1) + 1L]] <- group(paste(cc, arm, sep = " | "),
                                   x$condition == cc & d$randomization_group == arm, cc)
  hi <- d$randomization_group == "Intervention" & d$AI_HIGH == 1
  completed <- d$outcome_tte == 1 & d$tte_trigger == "Trial-triggered"
  completed <- completed %in% TRUE
  gs2 <- list(group("Completed", hi & completed), group("Did not complete", hi & !completed))
  gp <- lapply(seq_along(paths), function(i) group(paths[i], d$PATHWAY_INDEX == i))
  gv <- lapply(seq_along(paths), function(i) group(paths[i],
                                               d$PATHWAY_INDEX == i & d$outcome_vhd == 1))
  allgroups <- list(Table_1 = gmain, Supplemental_Table_1 = gs1,
                   Supplemental_Table_2 = gs2, Supplemental_Table_3 = gp,
                   Supplemental_Table_4 = gv)
  out <- list()
  for (nm in names(spec)) {
    s <- spec[[nm]]; groups <- allgroups[[nm]]
    n <- vapply(groups, function(g) sum(g$mask), integer(1))
    headers <- c("Characteristic", vapply(seq_along(groups), function(i)
      paste0(groups[[i]]$label, " (n = ", format(n[i], big.mark = ",", trim = TRUE), ")"), character(1)))
    has_tests <- nm %in% c("Supplemental_Table_2", "Supplemental_Table_3", "Supplemental_Table_4")
    if (has_tests) headers <- c(headers, "P value")
    cells <- matrix("", nrow(s), length(headers)); cells[, 1] <- s$label
    for (r in seq_len(nrow(s))) {
      key <- s$key[r]
      if (nzchar(key)) {
        for (j in seq_along(groups)) {
          # Preserve the intentional blank cells in Table 1.
          if (nm == "Table_1" && ((key %in% c("high", "low") && j > 2L) ||
              (key == "days_high" && !j %in% 3:4))) next
          cells[r, j + 1L] <- valve_table_summary(
            x[[key]][groups[[j]]$mask], s$kind[r], s$category[r])
        }
      }
      if (has_tests && nzchar(s$test_key[r])) {
        method <- if (s$test_kind[r] == "continuous") {
          if (nm == "Supplemental_Table_2") "student" else "anova"
        } else "fisher"
        z <- valve_table_test(x, groups, s$test_key[r], method, config)
        cells[r, ncol(cells)] <- valve_table_p(z$P_value, z$Chi_square_fallback)
      }
    }
    colnames(cells) <- headers
    out[[nm]] <- as.data.frame(cells, stringsAsFactors = FALSE, check.names = FALSE)
  }
  out
}

# Labels and row specifications for the five manuscript tables.
valve_table_schema <- function() {
  out <- list()
  out[["Table_1"]] <- data.frame(
    original_row = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69),
    label = c("ECG profile", "Stratification by AI-ECG", "High risk of VHD, n (%)", "Low risk of VHD, n (%)", "Days from most recent ECG, median (IQR)", "Days from most recent high-risk ECG, median (IQR)", "Demographics and prespecified clinical condition", "  Male, n (%)", "  Age (years, mean ± SD)", "Age group", "60-74 years, n (%)", "75-85 years, n (%)", "  Prespecified clinical condition subgroup", "High-priority echocardiography condition, n (%)", "Potential echocardiography condition, n (%)", "No prespecified echocardiography condition, n (%)", "Recent symptom (≤3 months)", "  Dyspnea, n (%)", "  Palpitations, n (%)", "  Angina or chest pain, n (%)", "  Syncope, n (%)", "  Orthopnea, n (%)", "  Edema, n (%)", "Clinical history and cardiovascular risk profile", "Aortic aneurysm, n (%)", "Marfan syndrome, n (%)", "  Atrial fibrillation", "    Recent-onset (≤3 years), n (%)", "Long-standing (>3 years), n (%)", "Heart failure", "Recent-onset (≤3 years), n (%)", "Long-standing (>3 years), n (%)", "Chronic kidney disease", "Recent-onset (≤3 years), n (%)", "Long-standing (>3 years), n (%)", "Coronary artery disease", "Recent-onset (≤3 years), n (%)", "Long-standing (>3 years), n (%)", "Myocardial infarction", "Recent-onset (≤3 years), n (%)", "Long-standing (>3 years), n (%)", "Prior left ventricular ejection fraction", "Severely reduced (≤40%), n (%)", "Mildly reduced (41-49%), n (%)", "Recent-onset markedly elevated NT-proBNP, n (%)†", "Moderate elevated NT-proBNP, n (%)‡", "Diabetes mellitus, n (%)", "Hypertension, n (%)", "Thyroid disorder, n (%)", "Permanent pacemaker, n (%)", "Elevated pulmonary artery systolic pressure, n (%)", "Marked septal hypertrophy, n (%)", "Large pericardial effusion, n (%)", "Other outcome-related histories", "Other structural or pericardial cardiac abnormalities, n (%)", "Recent psychiatric consultation (≤6 months), n (%)", "ACE inhibitor or ARB, n (%)", "Angiotensin receptor-neprilysin inhibitors, n (%)", "Beta blocker, n (%)", "Mineralocorticoid receptor antagonist, n (%)", "Sodium-glucose cotransporter 2 inhibitors, n (%)", "Loop diuretic, n (%)", "Calcium channel blocker, n (%)", "Thiazides, n (%)", "Alpha blocker, n (%)", "Other antihypertensive agents, n (%)", "Regular cardiology follow-up, n (%)", "Regular cardiovascular surgery follow-up, n (%)"),
    key = c("", "", "high", "low", "days_ecg", "days_high", "", "male", "age", "", "age60", "age75", "", "condition1", "condition2", "condition3", "", "dyspnea", "palpitations", "angina", "syncope", "orthopnea", "edema", "", "aneurysm", "marfan", "", "af_recent", "af_long", "", "hf_recent", "hf_long", "", "ckd_recent", "ckd_long", "", "cad_recent", "cad_long", "", "mi_recent", "mi_long", "", "ef40", "ef41", "nt_high", "nt_moderate", "dm", "htn", "thyroid", "pacemaker", "pasp", "septum", "effusion", "", "other", "psych", "ace", "arni", "bb", "mra", "sglt2", "loop", "ccb", "thiazide", "alpha", "otherhtn", "cv", "cvs"),
    kind = c("", "", "binary", "binary", "median", "median", "", "binary", "mean", "", "binary", "binary", "", "binary", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary"),
    category = c("", "", "1", "1", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "1", "", "1", "1", "1", "1", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"),
    test_key = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""),
    test_kind = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""),
    stringsAsFactors = FALSE, check.names = FALSE)
  out[["Supplemental_Table_1"]] <- data.frame(
    original_row = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63),
    label = c("ECG profile", "  Stratification by AI-ECG", "    High risk of VHD, n (%)", "    Low risk of VHD, n (%)", "  Days from most recent ECG, median (IQR)", "  Days from most recent high-risk ECG, median (IQR)", "Demographics", "  Male, n (%)", "  Age (years, mean ± SD)", "Age group", "60-74 years, n (%)", "75-85 years, n (%)", "High-priority echocardiography condition (≥1 following)", "  Aortic aneurysm, n (%)", "  Marfan syndrome, n (%)", "  Severely reduced LVEF (≤40%), n (%)", "  Recent-onset atrial fibrillation (≤3 years), n (%)", "  Recent-onset heart failure (≤3 years), n (%)", "  Recent-onset chronic kidney disease (≤3 years), n (%)", "  Recent-onset coronary artery disease (≤3 years), n (%)", "  Recent-onset myocardial infarction (≤3 years), n (%)", "  Recent dyspnea (≤3 months), n (%)", "  Recent palpitations (≤3 months), n (%)", "  Recent angina or chest pain (≤3 months), n (%)", "  Recent syncope (≤3 months), n (%)", "  Recent orthopnea (≤3 months), n (%)", "  Recent edema (≤3 months), n (%)", "  Recent-onset markedly elevated NT-proBNP, n (%)†", "Potential echocardiography condition (≥1 following)", "  Diabetes mellitus combining age ≥75 years, n (%)", "  Hypertension combining age ≥75 years, n (%)", "  Chronic kidney disease combining age ≥75 years, n (%)", "  Long-standing atrial fibrillation (>3 years), n (%)", "  Long-standing heart failure (>3 years), n (%)", "  Long-standing coronary artery disease (>3 years), n (%)", "  Long-standing myocardial infarction (>3 years), n (%)", "  Thyroid disorder, n (%)", "  Permanent pacemaker, n (%)", "  Moderate elevated NT-proBNP, n (%)‡", "  Mildly reduced LVEF (41-49%), n (%)", "  Elevated pulmonary artery systolic pressure, n (%)", "  Marked septal hypertrophy, n (%)", "  Large pericardial effusion, n (%)", "Other clinical history", "  Diabetes mellitus, n (%)", "  Hypertension, n (%)", "  Chronic kidney disease, n (%)", "Other outcome-related histories", "  Other structural or pericardial cardiac abnormalities, n (%)", "  Recent psychiatric consultation (≤6 months), n (%)", "  ACE inhibitor or ARB, n (%)", "  Angiotensin receptor-neprilysin inhibitors, n (%)", "  Beta blocker, n (%)", "  Mineralocorticoid receptor antagonist, n (%)", "  Sodium-glucose cotransporter 2 inhibitors, n (%)", "  Loop diuretic, n (%)", "  Calcium channel blocker, n (%)", "  Thiazides, n (%)", "  Alpha blocker, n (%)", "  Other antihypertensive agents, n (%)", "  Regular cardiology follow-up, n (%)", "  Regular cardiovascular surgery follow-up, n (%)"),
    key = c("", "", "high", "low", "days_ecg", "days_high", "", "male", "age", "", "age60", "age75", "", "aneurysm", "marfan", "ef40", "af_recent", "hf_recent", "ckd_recent", "cad_recent", "mi_recent", "dyspnea", "palpitations", "angina", "syncope", "orthopnea", "edema", "nt_high", "", "dm75", "htn75", "ckd75", "af_long", "hf_long", "cad_long", "mi_long", "thyroid", "pacemaker", "nt_moderate", "ef41", "pasp", "septum", "effusion", "", "dm", "htn", "ckd", "", "other", "psych", "ace", "arni", "bb", "mra", "sglt2", "loop", "ccb", "thiazide", "alpha", "otherhtn", "cv", "cvs"),
    kind = c("", "", "binary", "binary", "median", "median", "", "binary", "mean", "", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "", "binary", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary"),
    category = c("", "", "1", "1", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "", "1", "1", "1", "", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"),
    test_key = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""),
    test_kind = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""),
    stringsAsFactors = FALSE, check.names = FALSE)
  out[["Supplemental_Table_2"]] <- data.frame(
    original_row = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65),
    label = c("ECG profile", "  Days from most recent ECG, median (IQR)", "  Days from most recent high-risk ECG, median (IQR)", "Demographics and prespecified clinical condition", "  Male, n (%)", "  Age (years, mean ± SD)", "Age group", "60-74 years, n (%)", "75-85 years, n (%)", "  Prespecified clinical condition subgroup", "High-priority echocardiography condition, n (%)", "Potential echocardiography condition, n (%)", "No prespecified echocardiography condition, n (%)", "Recent symptom (≤3 months)", "  Dyspnea, n (%)", "  Palpitations, n (%)", "  Angina or chest pain, n (%)", "  Syncope, n (%)", "  Orthopnea, n (%)", "  Edema, n (%)", "Clinical history and cardiovascular risk profile", "  Aortic aneurysm, n (%)", "  Marfan syndrome, n (%)", "  Atrial fibrillation", "    Recent-onset (≤3 years), n (%)", "Long-standing (>3 years), n (%)", "  Heart failure", "Recent-onset (≤3 years), n (%)", "Long-standing (>3 years), n (%)", "  Chronic kidney disease", "Recent-onset (≤3 years), n (%)", "Long-standing (>3 years), n (%)", "  Coronary artery disease", "Recent-onset (≤3 years), n (%)", "Long-standing (>3 years), n (%)", "  Myocardial infarction", "Recent-onset (≤3 years), n (%)", "Long-standing (>3 years), n (%)", "  Prior left ventricular ejection fraction", "Severely reduced (≤40%), n (%)", "Mildly reduced (41-49%), n (%)", "  Recent-onset markedly elevated NT-proBNP, n (%)†", "  Moderate elevated NT-proBNP, n (%)‡", "  Diabetes mellitus, n (%)", "  Hypertension, n (%)", "  Thyroid disorder, n (%)", "  Permanent pacemaker, n (%)", "  Elevated pulmonary artery systolic pressure, n (%)", "  Marked septal hypertrophy, n (%)", "  Large pericardial effusion, n (%)", "Other outcome-related histories", "  Other structural or pericardial cardiac abnormalities, n (%)", "  Recent psychiatric consultation (≤6 months), n (%)", "  ACE inhibitor or ARB, n (%)", "  Angiotensin receptor-neprilysin inhibitors, n (%)", "  Beta blocker, n (%)", "  Mineralocorticoid receptor antagonist, n (%)", "  Sodium-glucose cotransporter 2 inhibitors, n (%)", "  Loop diuretic, n (%)", "  Calcium channel blocker, n (%)", "  Thiazides, n (%)", "  Alpha blocker, n (%)", "  Other antihypertensive agents, n (%)", "  Regular cardiology follow-up, n (%)", "  Regular cardiovascular surgery follow-up, n (%)"),
    key = c("", "days_ecg", "days_high", "", "male", "age", "", "age60", "age75", "", "condition1", "condition2", "condition3", "", "dyspnea", "palpitations", "angina", "syncope", "orthopnea", "edema", "", "aneurysm", "marfan", "", "af_recent", "af_long", "", "hf_recent", "hf_long", "", "ckd_recent", "ckd_long", "", "cad_recent", "cad_long", "", "mi_recent", "mi_long", "", "ef40", "ef41", "nt_high", "nt_moderate", "dm", "htn", "thyroid", "pacemaker", "pasp", "septum", "effusion", "", "other", "psych", "ace", "arni", "bb", "mra", "sglt2", "loop", "ccb", "thiazide", "alpha", "otherhtn", "cv", "cvs"),
    kind = c("", "median", "median", "", "binary", "mean", "", "binary", "binary", "", "binary", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary"),
    category = c("", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "1", "", "1", "1", "1", "1", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"),
    test_key = c("", "days_ecg", "days_high", "", "male", "age", "age_group", "", "", "condition", "", "", "", "", "dyspnea", "palpitations", "angina", "syncope", "orthopnea", "edema", "", "aneurysm", "marfan", "af_group", "", "", "hf_group", "", "", "ckd_group", "", "", "cad_group", "", "", "mi_group", "", "", "ef_group", "", "", "nt_high", "nt_moderate", "dm", "htn", "thyroid", "pacemaker", "pasp", "septum", "effusion", "", "other", "psych", "ace", "arni", "bb", "mra", "sglt2", "loop", "ccb", "thiazide", "alpha", "otherhtn", "cv", "cvs"),
    test_kind = c("", "continuous", "continuous", "", "fisher", "continuous", "fisher", "", "", "fisher", "", "", "", "", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "", "fisher", "fisher", "fisher", "", "", "fisher", "", "", "fisher", "", "", "fisher", "", "", "fisher", "", "", "fisher", "", "", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher"),
    stringsAsFactors = FALSE, check.names = FALSE)
  out[["Supplemental_Table_3"]] <- data.frame(
    original_row = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62),
    label = c("Demographics and prespecified clinical condition", "  Male, n (%)", "  Age (years, mean ± SD)", "Age group", "60-74 years, n (%)", "75-85 years, n (%)", "  Prespecified clinical condition subgroup", "High-priority echocardiography condition, n (%)", "Potential echocardiography condition, n (%)", "No prespecified echocardiography condition, n (%)", "Recent symptom (≤3 months)", "  Dyspnea, n (%)", "  Palpitations, n (%)", "  Angina or chest pain, n (%)", "  Syncope, n (%)", "  Orthopnea, n (%)", "  Edema, n (%)", "Clinical history and cardiovascular risk profile", "  Aortic aneurysm, n (%)", "  Marfan syndrome, n (%)", "  Atrial fibrillation", "    Recent-onset (≤3 years), n (%)", "    Long-standing (>3 years), n (%)", "  Heart failure", "    Recent-onset (≤3 years), n (%)", "    Long-standing (>3 years), n (%)", "  Chronic kidney disease", "    Recent-onset (≤3 years), n (%)", "    Long-standing (>3 years), n (%)", "  Coronary artery disease", "    Recent-onset (≤3 years), n (%)", "    Long-standing (>3 years), n (%)", "  Myocardial infarction", "    Recent-onset (≤3 years), n (%)", "    Long-standing (>3 years), n (%)", "  Prior left ventricular ejection fraction", "    Severely reduced (≤40%), n (%)", "    Mildly reduced (41-49%), n (%)", "  Recent-onset markedly elevated NT-proBNP, n (%)†", "  Moderate elevated NT-proBNP, n (%)‡", "  Diabetes mellitus, n (%)", "  Hypertension, n (%)", "  Thyroid disorder, n (%)", "  Permanent pacemaker, n (%)", "  Elevated pulmonary artery systolic pressure, n (%)", "  Marked septal hypertrophy, n (%)", "  Large pericardial effusion, n (%)", "Other outcome-related histories", "  Other structural or pericardial cardiac abnormalities, n (%)", "  Recent psychiatric consultation (≤6 months), n (%)", "  ACE inhibitor or ARB, n (%)", "  Angiotensin receptor-neprilysin inhibitors, n (%)", "  Beta blocker, n (%)", "  Mineralocorticoid receptor antagonist, n (%)", "  Sodium-glucose cotransporter 2 inhibitors, n (%)", "  Loop diuretic, n (%)", "  Calcium channel blocker, n (%)", "  Thiazides, n (%)", "  Alpha blocker, n (%)", "  Other antihypertensive agents, n (%)", "  Regular cardiology follow-up, n (%)", "  Regular cardiovascular surgery follow-up, n (%)"),
    key = c("", "male", "age", "", "age60", "age75", "", "condition1", "condition2", "condition3", "", "dyspnea", "palpitations", "angina", "syncope", "orthopnea", "edema", "", "aneurysm", "marfan", "", "af_recent", "af_long", "", "hf_recent", "hf_long", "", "ckd_recent", "ckd_long", "", "cad_recent", "cad_long", "", "mi_recent", "mi_long", "", "ef40", "ef41", "nt_high", "nt_moderate", "dm", "htn", "thyroid", "pacemaker", "pasp", "septum", "effusion", "", "other", "psych", "ace", "arni", "bb", "mra", "sglt2", "loop", "ccb", "thiazide", "alpha", "otherhtn", "cv", "cvs"),
    kind = c("", "binary", "mean", "", "binary", "binary", "", "binary", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary", "binary"),
    category = c("", "1", "1", "", "1", "1", "", "1", "1", "1", "", "1", "1", "1", "1", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"),
    test_key = c("", "male", "age", "age_group", "", "", "condition", "", "", "", "", "dyspnea", "palpitations", "angina", "syncope", "orthopnea", "edema", "", "aneurysm", "marfan", "af_group", "", "", "hf_group", "", "", "ckd_group", "", "", "cad_group", "", "", "mi_group", "", "", "ef_group", "", "", "nt_high", "nt_moderate", "dm", "htn", "thyroid", "pacemaker", "pasp", "septum", "effusion", "", "other", "psych", "ace", "arni", "bb", "mra", "sglt2", "loop", "ccb", "thiazide", "alpha", "otherhtn", "cv", "cvs"),
    test_kind = c("", "fisher", "continuous", "fisher", "", "", "fisher", "", "", "", "", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "", "fisher", "fisher", "fisher", "", "", "fisher", "", "", "fisher", "", "", "fisher", "", "", "fisher", "", "", "fisher", "", "", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher"),
    stringsAsFactors = FALSE, check.names = FALSE)
  out[["Supplemental_Table_4"]] <- data.frame(
    original_row = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56),
    label = c("Demographics and prespecified clinical condition", "  Male, n (%)", "  Age (years, mean ± SD)", "Age group", "60-74 years, n (%)", "75-85 years, n (%)", "  Prespecified clinical condition subgroup", "High-priority echocardiography condition, n (%)", "Potential echocardiography condition, n (%)", "No prespecified echocardiography condition, n (%)", "Echocardiographic report", "  Number of moderate/severe valvular heart disease", "1, n (%)", "≥2, n (%)", "  Aortic stenosis severity", "Minimal/Mild, n (%)", "Moderate, n (%)", "Severe, n (%)", "  Aortic regurgitation severity", "Minimal/Mild, n (%)", "Moderate, n (%)", "Severe, n (%)", "  Mitral regurgitation severity", "Minimal/Mild, n (%)", "Moderate, n (%)", "Severe, n (%)", "  Tricuspid regurgitation severity", "Minimal/Mild, n (%)", "Moderate, n (%)", "Severe, n (%)", "  Left ventricular ejection fraction (%)", "Left ventricular ejection fraction group", ">50%, n (%)", "41-50%, n (%)", "≤40%, n (%)", "  Left ventricular end-diastolic diameter (mm)", "Left ventricular end-diastolic diameter group ", "<60 mm, n (%)", "≥60 mm, n (%)", " Pulmonary artery systolic pressure (mmHg)", "Pulmonary artery systolic pressure group ", "<50 mmHg, n (%)", "≥50 mmHg, n (%)", "Pre-randomized symptom", "  Dyspnea, n (%)", "  Palpitations, n (%)", "  Angina or chest pain, n (%)", "  Syncope, n (%)", "  Orthopnea, n (%)", "Symptoms after VHD diagnosis and within 90 days after index time", "  Dyspnea, n (%)", "  Palpitations, n (%)", "  Angina or chest pain, n (%)", "  Syncope, n (%)", "  Orthopnea, n (%)", "  Edema, n (%)"),
    key = c("", "male", "age", "", "age_group", "age_group", "", "condition", "condition", "condition", "", "", "vhd_valve_count_category", "vhd_valve_count_category", "", "vhd_aortic_stenosis_severity", "vhd_aortic_stenosis_severity", "vhd_aortic_stenosis_severity", "", "vhd_aortic_regurgitation_severity", "vhd_aortic_regurgitation_severity", "vhd_aortic_regurgitation_severity", "", "vhd_mitral_regurgitation_severity", "vhd_mitral_regurgitation_severity", "vhd_mitral_regurgitation_severity", "", "vhd_tricuspid_regurgitation_severity", "vhd_tricuspid_regurgitation_severity", "vhd_tricuspid_regurgitation_severity", "vhd_lvef_percent", "", "vhd_lvef_category", "vhd_lvef_category", "vhd_lvef_category", "vhd_lvedd_mm", "", "vhd_lvedd_category", "vhd_lvedd_category", "vhd_pasp_mmhg", "", "vhd_pasp_category", "vhd_pasp_category", "", "dyspnea", "palpitations", "angina", "syncope", "orthopnea", "", "vhd_postindex_dyspnea", "vhd_postindex_palpitations", "vhd_postindex_angina_chest_pain", "vhd_postindex_syncope", "vhd_postindex_orthopnea", "vhd_postindex_edema"),
    kind = c("", "binary", "mean", "", "binary", "binary", "", "binary", "binary", "binary", "", "", "binary", "binary", "", "binary", "binary", "binary", "", "binary", "binary", "binary", "", "binary", "binary", "binary", "", "binary", "binary", "binary", "mean", "", "binary", "binary", "binary", "mean", "", "binary", "binary", "mean", "", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "", "binary", "binary", "binary", "binary", "binary", "binary"),
    category = c("", "1", "", "", "60-74", "75-85", "", "High-priority echocardiography condition", "Potential echocardiography condition", "No prespecified echocardiography condition", "", "", "1", ">=2", "", "Mininal/Mild", "Moderate", "Severe", "", "Mininal/Mild", "Moderate", "Severe", "", "Mininal/Mild", "Moderate", "Severe", "", "Mininal/Mild", "Moderate", "Severe", "", "", ">50%", "41-50%", "<=40%", "", "", "<60 mm", ">=60 mm", "", "", "<50 mmHg", ">=50 mmHg", "", "1", "1", "1", "1", "1", "", "1", "1", "1", "1", "1", "1"),
    test_key = c("", "male", "age", "age_group", "", "", "condition", "", "", "", "", "vhd_valve_count_category", "", "", "vhd_aortic_stenosis_severity", "", "", "", "vhd_aortic_regurgitation_severity", "", "", "", "vhd_mitral_regurgitation_severity", "", "", "", "vhd_tricuspid_regurgitation_severity", "", "", "", "vhd_lvef_percent", "vhd_lvef_category", "", "", "", "vhd_lvedd_mm", "vhd_lvedd_category", "", "", "vhd_pasp_mmhg", "vhd_pasp_category", "", "", "", "dyspnea", "palpitations", "angina", "syncope", "orthopnea", "", "vhd_postindex_dyspnea", "vhd_postindex_palpitations", "vhd_postindex_angina_chest_pain", "vhd_postindex_syncope", "vhd_postindex_orthopnea", "vhd_postindex_edema"),
    test_kind = c("", "fisher", "continuous", "fisher", "", "", "fisher", "", "", "", "", "fisher", "", "", "fisher", "", "", "", "fisher", "", "", "", "fisher", "", "", "", "fisher", "", "", "", "continuous", "fisher", "", "", "", "continuous", "fisher", "", "", "continuous", "fisher", "", "", "", "fisher", "fisher", "fisher", "fisher", "fisher", "", "fisher", "fisher", "fisher", "fisher", "fisher", "fisher"),
    stringsAsFactors = FALSE, check.names = FALSE)
  out
}

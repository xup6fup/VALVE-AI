
VALVE_SUBGROUPS <- list(
  "prespecified_clinical_triggers_for_tte" = list(name = "Prespecified clinical triggers for TTE", levels = c("Strong TTE triggers", "Potential TTE triggers", "No identified TTE triggers")),
  "subgroup_sex" = list(name = "Sex", levels = c("Female", "Male")),
  "subgroup_age" = list(name = "Age", levels = c("60-74 years", "75-85 years")),
  "subgroup_recent_relevant_symptoms" = list(name = "Recent relevant symptoms", levels = c("Yes", "No")),
  "subgroup_outpatient_cardiology_care" = list(name = "Outpatient cardiology care", levels = c("Yes", "No")),
  "subgroup_cardiovascular_prescription" = list(name = "History of cardiovascular medication prescriptions", levels = c("Yes", "No")),
  "subgroup_diabetes_mellitus_hypertension_histories" = list(name = "Histories of diabetes mellitus and hypertension", levels = c("None", "Diabetes mellitus only", "Hypertension only", "Both")),
  "subgroup_chronic_kidney_disease_history" = list(name = "History of chronic kidney disease", levels = c("Yes", "No")),
  "subgroup_cardiovascular_disease_history" = list(name = "History of cardiovascular diseases", levels = c("Yes", "No")),
  "subgroup_most_recent_ecg_ai_risk" = list(name = "AI risk status of the most recent ECG", levels = c("AI high-risk", "AI low-risk")),
  "subgroup_time_since_most_recent_ai_high_risk_ecg" = list(name = "Time since the most recent AI high-risk ECG", levels = c("ECG time ≤1 month", "ECG time >1 month")),
  "subgroup_ai_ecg_risk_profile_preceding_year" = list(name = "AI-ECG risk profile in the preceding year", levels = c("Only 1 ECG, high risk", "≥2 ECGs, ≤50% high risk", "≥2 ECGs, >50% high risk"))
)

VALVE_OUTCOMES <- list(
  "outcome_vhd" = list(label = "Moderate or severe VHD", section = "Primary outcome"),
  "outcome_aortic_stenosis" = list(label = "Aortic stenosis", section = "Secondary outcomes"),
  "outcome_aortic_regurgitation" = list(label = "Aortic regurgitation", section = "Secondary outcomes"),
  "outcome_mitral_regurgitation" = list(label = "Mitral regurgitation", section = "Secondary outcomes"),
  "outcome_tricuspid_regurgitation" = list(label = "Tricuspid regurgitation", section = "Secondary outcomes"),
  "outcome_other_structural_abnormality" = list(label = "Other structural abnormalities", section = "Secondary outcomes"),
  "outcome_new_psychiatric_consultation" = list(label = "New psychiatric consultation", section = "Safety outcome"),
  "outcome_tte_after_randomization" = list(label = "TTE performed after randomization", section = "TTE utilization"),
  "outcome_vhd_with_relevant_symptoms" = list(label = "VHD with relevant symptoms", section = "Other composite VHD outcomes"),
  "outcome_vhd_without_relevant_symptoms" = list(label = "VHD without relevant symptoms", section = "Other composite VHD outcomes"),
  "outcome_vhd_ntprobnp_ge900" = list(label = "VHD with elevated NT-proBNP", section = "Other composite VHD outcomes"),
  "outcome_left_sided_vhd" = list(label = "Left-sided VHD", section = "Other composite VHD outcomes"),
  "outcome_left_sided_vhd_lvedd_ge60" = list(label = "Left-sided VHD with enlarged LVEDD", section = "Other composite VHD outcomes"),
  "outcome_severe_vhd" = list(label = "Severe VHD", section = "Other composite VHD outcomes"),
  "outcome_severe_vhd_surgical_criteria" = list(label = "Severe VHD meeting surgical criteria", section = "Other composite VHD outcomes"),
  "outcome_vhd_composite_examinations" = list(label = "Composite examinations after VHD diagnosed", section = "Examinations after VHD diagnosed"),
  "outcome_vhd_chest_xray" = list(label = "Chest radiography", section = "Examinations after VHD diagnosed"),
  "outcome_vhd_ntprobnp_test" = list(label = "NT-proBNP testing", section = "Examinations after VHD diagnosed"),
  "outcome_vhd_high_sensitivity_troponin_i_test" = list(label = "High-sensitivity troponin I testing", section = "Examinations after VHD diagnosed"),
  "outcome_vhd_exercise_test" = list(label = "Exercise testing", section = "Examinations after VHD diagnosed"),
  "outcome_vhd_transesophageal_echo" = list(label = "Transesophageal echocardiography", section = "Examinations after VHD diagnosed"),
  "outcome_vhd_composite_managements" = list(label = "Composite managements after VHD diagnosed", section = "Managements after VHD diagnosed"),
  "outcome_vhd_cardiologist_referral" = list(label = "Cardiologist referral", section = "Managements after VHD diagnosed"),
  "outcome_vhd_pharmacotherapy_modification" = list(label = "Pharmacotherapy modification", section = "Managements after VHD diagnosed"),
  "outcome_vhd_hf_pharmacotherapy_modification" = list(label = "Heart failure pharmacotherapy modification", section = "Managements after VHD diagnosed"),
  "outcome_vhd_antihypertensive_pharmacotherapy_modification" = list(label = "Antihypertensive pharmacotherapy modification", section = "Managements after VHD diagnosed"),
  "outcome_vhd_cardiovascular_surgeon_referral" = list(label = "Cardiovascular surgeon referral", section = "Managements after VHD diagnosed"),
  "outcome_vhd_valve_intervention" = list(label = "Surgical or transcatheter valve intervention", section = "Managements after VHD diagnosed")
)

VALVE_PATHWAYS <- c("AI-ECG high risk | Intervention: trial-triggered", "AI-ECG high risk | Intervention: physician-triggered", "AI-ECG high risk | Control: physician-triggered", "AI-ECG low risk | Intervention: physician-triggered", "AI-ECG low risk | Control: physician-triggered")
for (.valve_key in names(VALVE_SUBGROUPS)) {
  VALVE_SUBGROUPS[[.valve_key]]$display_levels <-
    gsub(">=", "\u2265", gsub("<=", "\u2264",
      VALVE_SUBGROUPS[[.valve_key]]$levels, fixed = TRUE), fixed = TRUE)
}
rm(.valve_key)

valve_config <- function(...) {
  defaults <- list(
    encoding = "UTF-8-BOM",
    font_family = "sans",
    fisher_workspace = 2e8,
    table_fisher_workspace = 200000L
  )
  given <- list(...)
  unknown <- setdiff(names(given), names(defaults))
  if (length(unknown)) stop("Unknown configuration: ", paste(unknown, collapse = ", "))
  defaults[names(given)] <- given
  defaults
}

valve_read_input <- function(input, config = valve_config()) {
  d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE,
                na.strings = c("", "NA", "NaN", "NULL"),
                fileEncoding = config$encoding)
  if (!nrow(d)) stop("The analysis requires participant records, not the header-only template.")
  if (!"case_id" %in% names(d) || anyNA(d$case_id) ||
      any(!nzchar(trimws(d$case_id))) || anyDuplicated(d$case_id))
    stop("case_id must be present, nonmissing, and unique.")
  if (anyNA(d$randomization_group) ||
      any(!d$randomization_group %in% c("Control", "Intervention")))
    stop("randomization_group must contain Control or Intervention.")
  d$randomization_group <- factor(d$randomization_group,
                                  levels = c("Control", "Intervention"))

  # Convenience variables for the analysis; the supplied clinical and outcome
  # indicators are not reclassified or reconstructed.
  d$AI_HIGH <- d$ai_ecg_risk=="High"; d$AI_RISK <- d$ai_ecg_risk
  d$CONDITION <- d$prespecified_clinical_triggers_for_tte
  d$DAYS_RECENT_ECG <- d$days_since_last_ecg
  d$DAYS_RECENT_HIGH_RISK_ECG <- ifelse(d$AI_HIGH,d$days_since_last_high_risk_ecg,NA_real_)
  ecg_lists <- lapply(as.character(d$ecg_days_before_index),function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return(numeric())
    values <- suppressWarnings(as.numeric(trimws(strsplit(s,";",fixed=TRUE)[[1]])))
    values
  })
  d$n_ECG <- vapply(ecg_lists,function(v)if(length(v))length(v) else NA_integer_,0L)
  d$TTE_SOURCE <- ifelse(is.na(d$tte_trigger),"No TTE",d$tte_trigger)
  d$OUTCOME_TTE_TRIAL <- as.integer(d$outcome_tte_after_randomization==1 & d$TTE_SOURCE=="Trial-triggered")
  d$OUTCOME_TTE_PHYSICIAN <- as.integer(d$outcome_tte_after_randomization==1 & d$TTE_SOURCE=="Physician-triggered")
  d$PATHWAY_INDEX <- NA_integer_
  is_tte <- d$outcome_tte_after_randomization %in% 1
  d$PATHWAY_INDEX[is_tte & d$AI_HIGH & d$randomization_group=="Intervention" & d$TTE_SOURCE=="Trial-triggered"] <- 1L
  d$PATHWAY_INDEX[is_tte & d$AI_HIGH & d$randomization_group=="Intervention" & d$TTE_SOURCE=="Physician-triggered"] <- 2L
  d$PATHWAY_INDEX[is_tte & d$AI_HIGH & d$randomization_group=="Control" & d$TTE_SOURCE=="Physician-triggered"] <- 3L
  d$PATHWAY_INDEX[is_tte & !d$AI_HIGH & d$randomization_group=="Intervention" & d$TTE_SOURCE=="Physician-triggered"] <- 4L
  d$PATHWAY_INDEX[is_tte & !d$AI_HIGH & d$randomization_group=="Control" & d$TTE_SOURCE=="Physician-triggered"] <- 5L
  d$PATHWAY <- VALVE_PATHWAYS[d$PATHWAY_INDEX]
  for (i in 1:8) d[[paste0("SG",i)]] <- d[[c("subgroup_sex", "subgroup_age", "subgroup_recent_relevant_symptoms", "subgroup_outpatient_cardiology_care", "subgroup_cardiovascular_prescription", "subgroup_diabetes_mellitus_hypertension_histories", "subgroup_chronic_kidney_disease_history", "subgroup_cardiovascular_disease_history")[i]]]
  for (i in 1:3) d[[paste0("ECG_SG",i)]] <- d[[c("subgroup_most_recent_ecg_ai_risk", "subgroup_time_since_most_recent_ai_high_risk_ecg", "subgroup_ai_ecg_risk_profile_preceding_year")[i]]]
  attr(d, "ecg_lists") <- ecg_lists
  d
}

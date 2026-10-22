
# Age subgroup labels: 60-74 years and 75-85 years.
VALVE_SUBGROUPS <- list(
  "clinical_condition" = list(name = "Prespecified clinical condition", levels = c("High-priority echocardiography condition", "Potential echocardiography condition", "No prespecified echocardiography condition")),
  "subgroup_sex" = list(name = "Sex", levels = c("Female", "Male")),
  "subgroup_age" = list(name = "Age", levels = c("60-74 years", "75-85 years")),
  "subgroup_recent_symptoms" = list(name = "Recent qualifying symptoms", levels = c("Any qualifying symptom", "No qualifying symptom")),
  "subgroup_cardiology_followup" = list(name = "Regular cardiology follow-up", levels = c("Regular follow-up", "No regular follow-up")),
  "subgroup_cardiovascular_prescription" = list(name = "History of cardiovascular medication prescriptions", levels = c("Any qualifying recent prescription", "No qualifying recent prescription")),
  "subgroup_diabetes_hypertension" = list(name = "Diabetes mellitus and hypertension phenotype", levels = c("Neither diagnosis recorded", "Diabetes mellitus only", "Hypertension only", "Both diagnoses recorded")),
  "subgroup_ckd_history" = list(name = "Recorded chronic kidney disease history", levels = c("Recorded history", "No recorded history")),
  "subgroup_cardiac_history" = list(name = "Recorded cardiac disease", levels = c("Any recorded disease", "No recorded history")),
  "subgroup_latest_ecg_risk" = list(name = "Most recent ECG risk classification", levels = c("High-risk most recent ECG", "Low-risk most recent ECG")),
  "subgroup_high_risk_ecg_recency" = list(name = "Time from most recent high-risk ECG to index time", levels = c("<=1 month", ">1 month")),
  "subgroup_ecg_high_risk_pattern" = list(name = "ECG distribution in the preceding year", levels = c("Exactly 1 ECG, high risk", ">=2 ECGs, <=50% high risk", ">=2 ECGs, >50% high risk"))
)

# Post-VHD outcomes require both: after VHD diagnosis and within 90 days
# after index time. The supplied final indicators are used directly.
VALVE_OUTCOMES <- list(
  "outcome_vhd" = list(label = "Moderate or severe VHD", section = "Primary outcome"),
  "outcome_aortic_stenosis" = list(label = "Aortic stenosis", section = "Secondary outcomes"),
  "outcome_aortic_regurgitation" = list(label = "Aortic regurgitation", section = "Secondary outcomes"),
  "outcome_mitral_regurgitation" = list(label = "Mitral regurgitation", section = "Secondary outcomes"),
  "outcome_tricuspid_regurgitation" = list(label = "Tricuspid regurgitation", section = "Secondary outcomes"),
  "outcome_other_structural_abnormality" = list(label = "Other structural or pericardial abnormalities", section = "Secondary outcomes"),
  "outcome_new_psychiatry_visit" = list(label = "New psychiatric consultation", section = "Safety outcome"),
  "outcome_tte" = list(label = "TTE for any reason", section = "TTE utilization"),
  "outcome_vhd_with_symptoms" = list(label = "VHD with relevant symptoms", section = "Other composite VHD outcomes"),
  "outcome_vhd_without_recorded_symptoms" = list(label = "VHD without relevant symptoms", section = "Other composite VHD outcomes"),
  "outcome_vhd_ntprobnp_ge900" = list(label = "VHD with elevated NT-proBNP", section = "Other composite VHD outcomes"),
  "outcome_left_sided_vhd" = list(label = "Left-sided VHD", section = "Other composite VHD outcomes"),
  "outcome_left_sided_vhd_lvedd_ge60" = list(label = "Left-sided VHD with enlarged LVEDD", section = "Other composite VHD outcomes"),
  "outcome_severe_vhd" = list(label = "Severe VHD", section = "Other composite VHD outcomes"),
  "outcome_severe_vhd_surgical_criteria" = list(label = "Severe VHD meeting surgical criteria", section = "Other composite VHD outcomes"),
  "outcome_vhd_any_examination" = list(label = "Composite examination action", section = "Post-VHD examination outcomes"),
  "outcome_vhd_chest_xray" = list(label = "Chest radiography", section = "Post-VHD examination outcomes"),
  "outcome_vhd_ntprobnp_test" = list(label = "NT-proBNP testing", section = "Post-VHD examination outcomes"),
  "outcome_vhd_troponin_test" = list(label = "High-sensitivity troponin I testing", section = "Post-VHD examination outcomes"),
  "outcome_vhd_exercise_test" = list(label = "Exercise testing", section = "Post-VHD examination outcomes"),
  "outcome_vhd_transesophageal_echo" = list(label = "Transesophageal echocardiography", section = "Post-VHD examination outcomes"),
  "outcome_vhd_any_management" = list(label = "Composite management action", section = "Post-VHD management outcomes"),
  "outcome_vhd_cardiology_referral" = list(label = "Cardiology referral", section = "Post-VHD management outcomes"),
  "outcome_vhd_any_medication_change" = list(label = "Pharmacotherapy modification", section = "Post-VHD management outcomes"),
  "outcome_vhd_hf_medication_change" = list(label = "Heart failure pharmacotherapy modification", section = "Post-VHD management outcomes"),
  "outcome_vhd_antihypertensive_change" = list(label = "Antihypertensive pharmacotherapy modification", section = "Post-VHD management outcomes"),
  "outcome_vhd_surgery_referral" = list(label = "Cardiovascular surgery referral", section = "Post-VHD management outcomes"),
  "outcome_vhd_valve_intervention" = list(label = "Surgical or transcatheter valve intervention", section = "Post-VHD management outcomes")
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
  if (anyNA(d$randomization_group) ||
      any(!d$randomization_group %in% c("Control", "Intervention")))
    stop("randomization_group must contain Control or Intervention.")
  d$randomization_group <- factor(d$randomization_group,
                                  levels = c("Control", "Intervention"))

  # Convenience variables for the analysis; the supplied clinical and outcome
  # indicators are not reclassified or reconstructed.
  d$AI_HIGH <- d$ai_ecg_risk=="High"; d$AI_RISK <- d$ai_ecg_risk
  d$CONDITION <- d$clinical_condition
  d$DAYS_RECENT_ECG <- d$days_since_last_ecg
  d$DAYS_RECENT_HIGH_RISK_ECG <- ifelse(d$AI_HIGH,d$days_since_last_high_risk_ecg,NA_real_)
  ecg_lists <- lapply(as.character(d$ecg_days_before_index),function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return(numeric())
    values <- suppressWarnings(as.numeric(trimws(strsplit(s,";",fixed=TRUE)[[1]])))
    values
  })
  d$n_ECG <- vapply(ecg_lists,function(v)if(length(v))length(v) else NA_integer_,0L)
  d$TTE_SOURCE <- ifelse(is.na(d$tte_trigger),"No TTE",d$tte_trigger)
  d$OUTCOME_TTE_TRIAL <- as.integer(d$outcome_tte==1 & d$TTE_SOURCE=="Trial-triggered")
  d$OUTCOME_TTE_PHYSICIAN <- as.integer(d$outcome_tte==1 & d$TTE_SOURCE=="Physician-triggered")
  d$PATHWAY_INDEX <- NA_integer_
  is_tte <- d$outcome_tte %in% 1
  d$PATHWAY_INDEX[is_tte & d$AI_HIGH & d$randomization_group=="Intervention" & d$TTE_SOURCE=="Trial-triggered"] <- 1L
  d$PATHWAY_INDEX[is_tte & d$AI_HIGH & d$randomization_group=="Intervention" & d$TTE_SOURCE=="Physician-triggered"] <- 2L
  d$PATHWAY_INDEX[is_tte & d$AI_HIGH & d$randomization_group=="Control" & d$TTE_SOURCE=="Physician-triggered"] <- 3L
  d$PATHWAY_INDEX[is_tte & !d$AI_HIGH & d$randomization_group=="Intervention" & d$TTE_SOURCE=="Physician-triggered"] <- 4L
  d$PATHWAY_INDEX[is_tte & !d$AI_HIGH & d$randomization_group=="Control" & d$TTE_SOURCE=="Physician-triggered"] <- 5L
  d$PATHWAY <- VALVE_PATHWAYS[d$PATHWAY_INDEX]
  for (i in 1:8) d[[paste0("SG",i)]] <- d[[c("subgroup_sex", "subgroup_age", "subgroup_recent_symptoms", "subgroup_cardiology_followup", "subgroup_cardiovascular_prescription", "subgroup_diabetes_hypertension", "subgroup_ckd_history", "subgroup_cardiac_history")[i]]]
  for (i in 1:3) d[[paste0("ECG_SG",i)]] <- d[[c("subgroup_latest_ecg_risk", "subgroup_high_risk_ecg_recency", "subgroup_ecg_high_risk_pattern")[i]]]
  attr(d, "ecg_lists") <- ecg_lists
  d
}

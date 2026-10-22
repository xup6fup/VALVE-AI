# VALVE-AI analysis project

This R project reproduces the tables and figures for **Repurposing prior electrocardiograms for AI-guided detection of valvular heart disease: a VALVE-AI randomized clinical trial**. The final analysis includes **7,016 participants and 114 variables**: 3,518 participants assigned to intervention and 3,498 assigned to control. Each row represents one randomized participant.

The analysis starts directly from the prepared de-identified dataset with readable variable names. The variable definitions are provided below. Running the project produces **22 individual PDF figures and five CSV tables**.

## Data access

**Participant-level data are not publicly distributed** because the combination of 114 variables carries a potential re-identification risk. The public repository provides the analysis code, variable definitions, and a header-only input template.

The de-identified analysis dataset may be requested from the corresponding author (Chin Lin, xup6fup0629@gmail.com), with access subject to approval following a simplified ethics review. Full individual-patient data that can be linked to longitudinal electronic health records, and raw ECG signals, require formal review and approval by the ethics committee of Tri-Service General Hospital. These full records and raw signals are not required to run this project.

## Run the analyses

Use R with Cairo PDF support. Install the required analysis package once:

```r
install.packages("logistf", repos = "https://cloud.r-project.org")
```

Place the approved 114-variable dataset at `data/VALVE_AI.csv`. Open `VALVE-AI.Rproj` in RStudio and run:

```r
source("run_all.R")
```

RStudio is optional. From the project directory, the equivalent command is:

```sh
Rscript run_all.R
```

All calculations and graphics are implemented in R. The analysis uses `logistf` for Firth regression and the standard R packages `stats`, `graphics`, `grDevices`, `grid`, and `utils` for the remaining calculations, plotting, and file input/output. No packages are installed automatically, and the analysis can run offline once its dependencies are installed. The header-only template contains no participant records and cannot reproduce the findings.

The only generated files are the **22 PDFs in `output/figures/`** and the **five CSVs in `output/tables/`**. Repeated runs overwrite these files.

## Project contents

| Path | Contents |
|---|---|
| `VALVE-AI.Rproj` | RStudio project |
| `run_all.R` | One-command analysis entry point |
| `data/README.md` | Data-access and input instructions |
| `data/VALVE_AI_template.csv` | Header-only template containing the 114 readable variable names |
| `R/01_data.R` | Input loading, analysis groups, and plotting categories |
| `R/02_statistics.R` | Logistic and Firth regression, interaction tests, Wilson intervals, and Fisher tests |
| `R/03_tables.R` | Calculations and row specifications for the five tables |
| `R/04_figures.R` | Figure calculations and plotting functions |
| `R/05_exports.R` | Export of the 22 PDFs and five CSVs |
| `output/figures/` | 22 individual PDF figures |
| `output/tables/` | Five CSV tables |

## Figure and table outputs

The project produces main Figures 2–5, Extended Data Figures 1–10, Supplemental Figures 1–8, Table 1, and Supplemental Tables 1–4.

| Filename | Contents |
|---|---|
| `Fig_02.pdf` | Primary, secondary and safety outcomes by randomized group |
| `Fig_03.pdf` | TTE utilization and VHD detection by AI-ECG risk |
| `Fig_04.pdf` | Other prespecified composite VHD outcomes |
| `Fig_05.pdf` | Post-VHD examination and management outcomes |
| `Ex_Data_Fig_01.pdf` | Timing and number of available ECGs |
| `Ex_Data_Fig_02.pdf` | Primary outcome across clinical and demographic subgroups, all participants |
| `Ex_Data_Fig_03.pdf` | Primary outcome across clinical and demographic subgroups, AI-ECG high risk |
| `Ex_Data_Fig_04.pdf` | TTE utilization and primary outcome across ECG subgroups |
| `Ex_Data_Fig_05.pdf` | TTE pathways and diagnostic yields |
| `Ex_Data_Fig_06.pdf` | Clinical-condition composition and diagnostic yield across five TTE pathways |
| `Ex_Data_Fig_07.pdf` | Trial-triggered TTE diagnostic yield across subgroups |
| `Ex_Data_Fig_08.pdf` | Other composite VHD outcomes across five TTE pathways |
| `Ex_Data_Fig_09.pdf` | Conditional examination proportions across five TTE pathways |
| `Ex_Data_Fig_10.pdf` | Conditional management proportions across five TTE pathways |
| `Supplemental_Fig_01.pdf` | Diagnostic yields by sex, age, symptoms and cardiology follow-up |
| `Supplemental_Fig_02.pdf` | Diagnostic yields by medications and comorbidities |
| `Supplemental_Fig_03.pdf` | Conditional examination composite by clinical condition, sex and age |
| `Supplemental_Fig_04.pdf` | Conditional examination composite by symptoms, cardiology follow-up and medications |
| `Supplemental_Fig_05.pdf` | Conditional examination composite by comorbidities |
| `Supplemental_Fig_06.pdf` | Conditional management composite by clinical condition, sex and age |
| `Supplemental_Fig_07.pdf` | Conditional management composite by symptoms, cardiology follow-up and medications |
| `Supplemental_Fig_08.pdf` | Conditional management composite by comorbidities |
| `Table_1.csv` | Baseline characteristics by randomized group and AI-ECG risk |
| `Supplemental_Table_1.csv` | Baseline characteristics by prespecified clinical condition |
| `Supplemental_Table_2.csv` | High-risk intervention participants by trial-triggered TTE completion |
| `Supplemental_Table_3.csv` | Baseline characteristics of TTE recipients across five pathways |
| `Supplemental_Table_4.csv` | Characteristics of detected VHD cases across five pathways |

## Statistical methods

All randomized participants are analyzed in their assigned groups, with control as the regression reference. Each outcome uses its nonmissing eligible observations as the denominator. The supplied final values are used directly: 1 denotes an event, 0 denotes no event, and blank values remain missing or not applicable. No imputation or participant deduplication is performed.

Treatment effects are estimated using unadjusted binomial logistic regression with Wald 95% confidence intervals and two-sided coefficient P values. If an estimable comparison in an outcome/subgroup block has zero events in one arm, the block uses Firth penalized logistic regression with profile penalized-likelihood confidence intervals and penalized likelihood-ratio P values. Comparisons with zero events in both arms are not estimable and do not, by themselves, switch other comparisons to Firth regression.

Treatment-by-subgroup interactions are evaluated by a global likelihood-ratio comparison of interaction and additive models. When Firth regression is required, the interaction coefficients are jointly tested within the full penalized model. An empty group-by-subgroup cell or a stratum with zero events in both arms yields a nonestimable interaction.

Diagnostic yields and conditional action proportions use two-sided Wilson score 95% confidence intervals without continuity correction. Pathway-yield comparisons use two-sided Fisher exact tests; multilevel comparisons use global Fisher–Freeman–Halton tests. An empty denominator is not estimable. When both groups have observed participants but no events, P = 1.

Table 1 and Supplemental Table 1 are descriptive. Supplemental Table 2 uses two-sided Student t tests with pooled variance for continuous variables, including ECG intervals displayed as median (IQR). Supplemental Tables 3 and 4 use one-way ANOVA for continuous variables. Each multilevel categorical variable receives one global test, shown on its parent row.

Categorical comparisons in Supplemental Tables 2–4 first use `stats::fisher.test(..., workspace = 200000, simulate.p.value = FALSE)`. Only a computation error triggers `stats::chisq.test(..., simulate.p.value = FALSE, correct = TRUE)`: Pearson's asymptotic chi-squared test, with Yates's continuity correction for 2 × 2 tables. These fallback P values are marked with §. Small expected counts do not, by themselves, trigger the fallback. Constant categorical variables retain a zero-count complement for the Fisher test, yielding P = 1.000. No simulated P values are used.

The chi-squared fallback applies only to the categorical comparisons in Supplemental Tables 2–4. Diagnostic-yield and conditional-action comparisons continue to use Fisher exact tests; a failed exact calculation is reported as not estimable.

### Table and figure notation

Table values are mean ± SD, median (IQR), or n (%), as indicated by the row labels, and use available observations for each variable. ECG interval medians and IQRs are displayed as integers. In the tables, † denotes NT-proBNP ≥900 pg/mL within three years, ‡ denotes NT-proBNP ≥300 pg/mL in sinus rhythm or ≥600 pg/mL in atrial fibrillation, and § marks a Pearson chi-squared fallback P value. Prior LVEF categories include ≤40%, 41–49%, and no documented reduction; the last category does not imply a measured LVEF ≥50%. In regression figures, † identifies Firth-based inference. NE denotes not estimable, and Ref. denotes the reference pathway.

SD, standard deviation; IQR, interquartile range; AI-ECG, artificial intelligence-enabled electrocardiogram; TTE, transthoracic echocardiography; VHD, valvular heart disease; LVEF, left ventricular ejection fraction; LVEDD, left ventricular end-diastolic diameter; PASP, pulmonary artery systolic pressure; NT-proBNP, N-terminal pro-B-type natriuretic peptide.

## Variable definitions

The following definitions correspond to the prepared input columns. Binary clinical flags use 1 for a recorded qualifying condition and 0 when no qualifying condition is recorded; 0 does not necessarily establish the clinical absence of disease. Blank values remain missing or not applicable. Column names must match the input template, and categorical values must match those listed below; the supplied `Mininal/Mild` code is displayed as `Minimal/Mild` in tables.

| Variable | Type | Coding | Definition |
|---|---|---|---|
| `randomization_group` | arm | Control  /  Intervention | Randomized allocation. Intervention denotes the AI-ECG-guided strategy; Control denotes usual care and is the regression reference. |
| `ai_ecg_risk` | categorical | High  /  Low | Precomputed baseline AI-ECG risk: any qualifying high-risk ECG in the preceding year classifies the participant as High. |
| `days_since_last_ecg` | numeric | Nonnegative days | Days from the most recent available ECG to the index time. |
| `days_since_last_high_risk_ecg` | numeric | Nonnegative days; blank = not applicable or unavailable | Days from the most recent high-risk ECG to the index time. Applicable to AI-ECG High participants. |
| `ecg_days_before_index` | ecg_day_list | Semicolon-separated nonnegative integer days; blank = unavailable | Semicolon-separated days from each available ECG to index time, for example 12;25;100. Multiple ECGs on the same day remain separate entries. Repeated day values are not removed. |
| `sex` | categorical | female  /  male | Recorded baseline sex category. |
| `age_years` | numeric | Numeric years | Age in years at the index time. |
| `age_group` | categorical | 60-74  /  75-85 | Prespecified age categories at index time: 60–74 years and 75–85 years, using the supplied category values. |
| `high_priority_tte_condition` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | At least one high-priority TTE condition, as indicated by the supplied clinical flag. |
| `history_aortic_aneurysm` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Documented history of aortic aneurysm. |
| `history_marfan_syndrome` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Documented history of Marfan syndrome. |
| `history_lvef_le40` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Prior left ventricular ejection fraction (LVEF) at most 40%. |
| `atrial_fibrillation_within_3y` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Atrial fibrillation documented within the preceding 3 years. |
| `heart_failure_within_3y` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Heart failure documented within the preceding 3 years. |
| `ckd_within_3y` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Chronic kidney disease documented within the preceding 3 years. |
| `coronary_disease_within_3y` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Coronary artery disease documented within the preceding 3 years. |
| `myocardial_infarction_within_3y` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Myocardial infarction documented within the preceding 3 years. |
| `ntprobnp_ge900_within_3y` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Elevated NT-proBNP within the preceding 3 years; the clinical threshold is at least 900 pg/mL. |
| `dyspnea_within_3m` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Dyspnea within the preceding 3 months. |
| `palpitations_within_3m` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Palpitations within the preceding 3 months. |
| `angina_chest_pain_within_3m` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Angina or chest pain within the preceding 3 months. |
| `syncope_within_3m` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Syncope within the preceding 3 months. |
| `orthopnea_within_3m` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Orthopnea within the preceding 3 months. |
| `edema_within_3m` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Edema within the preceding 3 months. |
| `potential_tte_condition` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | At least one potential TTE condition. This indicator may coexist with high_priority_tte_condition; the mutually exclusive hierarchy is applied only in clinical_condition. |
| `diabetes_age_ge75` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Diabetes mellitus with age at least 75 years. |
| `hypertension_age_ge75` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Hypertension with age at least 75 years. |
| `ckd_over_3y_age_ge75` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Chronic kidney disease documented more than 3 years previously, with age at least 75 years. |
| `history_thyroid_disorder` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Documented history of thyroid disease. |
| `atrial_fibrillation_over_3y` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Atrial fibrillation documented more than 3 years previously. |
| `heart_failure_over_3y` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Heart failure documented more than 3 years previously. |
| `coronary_disease_over_3y` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Coronary artery disease documented more than 3 years previously. |
| `myocardial_infarction_over_3y` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Myocardial infarction documented more than 3 years previously. |
| `history_pacemaker` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Documented history of a permanent pacemaker. |
| `history_elevated_ntprobnp` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Prior elevated NT-proBNP, using thresholds of at least 300 pg/mL in sinus rhythm or at least 600 pg/mL in atrial fibrillation. |
| `history_lvef_below50` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Prior reduced LVEF (<50%), represented by the supplied flag. |
| `history_elevated_pasp` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Prior elevated pulmonary artery systolic pressure (PASP>50 mmHg), represented by the supplied flag. |
| `history_lv_hypertrophy` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Prior left ventricular hypertrophy (interventricular septal thickness ≥15 mm), represented by the supplied flag. |
| `history_large_pericardial_effusion` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Prior large pericardial effusion (≥20 mm), represented by the supplied flag. |
| `history_diabetes` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Baseline history of diabetes mellitus without an age restriction. |
| `history_hypertension` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Baseline history of hypertension without an age restriction. |
| `ckd_over_3y` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Chronic kidney disease documented more than 3 years previously, without an age restriction. This is not equivalent to any prior history of chronic kidney disease. |
| `history_other_structural_abnormality` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Prior atrial septal defect, ventricular septal defect, cardiac tamponade or pericardial effusion. |
| `psychiatry_visit_within_180d` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Psychiatric contact within the 180 days before baseline; this defines an exclusion from the safety outcome risk set. |
| `baseline_acei_arb` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Baseline prescription history of an angiotensin-converting enzyme inhibitor or angiotensin receptor blocker. |
| `baseline_arni` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Baseline prescription history of an angiotensin receptor-neprilysin inhibitor. |
| `baseline_beta_blocker` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Baseline prescription history of a beta blocker. |
| `baseline_mra` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Baseline prescription history of a mineralocorticoid receptor antagonist. |
| `baseline_sglt2_inhibitor` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Baseline prescription history of a sodium-glucose cotransporter 2 inhibitor. |
| `baseline_loop_diuretic` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Baseline prescription history of a loop diuretic. |
| `baseline_calcium_channel_blocker` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Baseline prescription history of a calcium channel blocker. |
| `baseline_thiazide` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Baseline prescription history of a thiazide diuretic. |
| `baseline_alpha_blocker` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Baseline prescription history of an alpha blocker. |
| `baseline_other_antihypertensive` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Baseline prescription history of another antihypertensive medication. |
| `regular_cardiology_followup` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | At least two cardiology outpatient records from 365 days before the index time through the index time. |
| `regular_cardiovascular_surgery_followup` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | At least two cardiovascular surgery outpatient records from 365 days before the index time through the index time. |
| `clinical_condition` | categorical | High-priority echocardiography condition  /  Potential echocardiography condition  /  No prespecified echocardiography condition | Three supplied mutually exclusive clinical conditions: high-priority criteria take precedence, followed by potential criteria, with neither criterion forming the remaining group. Category labels are High-priority echocardiography condition, Potential echocardiography condition and No prespecified echocardiography condition. |
| `subgroup_sex` | categorical | Female  /  Male | Supplied sex subgroup. |
| `subgroup_age` | categorical | 60-74 years  /  75-85 years | Prespecified age subgroup at index time: 60–74 years versus 75–85 years. |
| `subgroup_recent_symptoms` | categorical | Any qualifying symptom  /  No qualifying symptom | Any qualifying dyspnea, palpitations, angina/chest pain, syncope, orthopnea or edema within the 3 months before the index time, according to the six supplied baseline symptom flags. |
| `subgroup_cardiology_followup` | categorical | Regular follow-up  /  No regular follow-up | Regular cardiology follow-up, based on at least two cardiology outpatient records in the year before the index time. |
| `subgroup_cardiovascular_prescription` | categorical | Any qualifying recent prescription  /  No qualifying recent prescription | Any qualifying baseline prescription of ACEi/ARB, ARNI, beta blocker, MRA, SGLT2 inhibitor, loop diuretic, calcium channel blocker, thiazide, alpha blocker or another antihypertensive. |
| `subgroup_diabetes_hypertension` | categorical | Neither diagnosis recorded  /  Diabetes mellitus only  /  Hypertension only  /  Both diagnoses recorded | Four mutually exclusive categories combining diabetes mellitus and hypertension without an age restriction. |
| `subgroup_ckd_history` | categorical | Recorded history  /  No recorded history | Any baseline chronic kidney disease history. Unlike ckd_over_3y, this is not restricted to disease documented more than 3 years previously. |
| `subgroup_cardiac_history` | categorical | Any recorded disease  /  No recorded history | Any baseline history of atrial fibrillation, heart failure or coronary artery disease. |
| `subgroup_latest_ecg_risk` | categorical | High-risk most recent ECG  /  Low-risk most recent ECG | Risk category of the most recent ECG among AI-ECG High participants. |
| `subgroup_high_risk_ecg_recency` | categorical | <=1 month  /  >1 month | Recency of the most recent high-risk ECG among AI-ECG High participants: at most 30 days versus more than 30 days. |
| `subgroup_ecg_high_risk_pattern` | categorical | Exactly 1 ECG, high risk  /  >=2 ECGs, <=50% high risk  /  >=2 ECGs, >50% high risk | Supplied ECG pattern in the preceding year among AI-ECG High participants: exactly one eligible ECG, classified as high risk; at least two ECGs with at most 50% high risk; or at least two ECGs with more than 50% high risk. Exactly 50% belongs to the <=50% category. |
| `outcome_vhd` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Newly detected moderate or severe aortic stenosis, aortic regurgitation, mitral regurgitation or tricuspid regurgitation within 90 days after randomization. |
| `outcome_aortic_stenosis` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Newly detected moderate or severe aortic stenosis within 90 days after randomization. |
| `outcome_aortic_regurgitation` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Newly detected moderate or severe aortic regurgitation within 90 days after randomization. |
| `outcome_mitral_regurgitation` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Newly detected moderate or severe mitral regurgitation within 90 days after randomization. |
| `outcome_tricuspid_regurgitation` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Newly detected moderate or severe tricuspid regurgitation within 90 days after randomization. |
| `outcome_other_structural_abnormality` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Newly detected atrial septal defect, ventricular septal defect, cardiac tamponade or large pericardial effusion (at least 20 mm) within 90 days. |
| `outcome_new_psychiatry_visit` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | New psychiatric consultation within 90 days. Participants with psychiatric contact in the 180 days before baseline have missing/not-applicable values and are excluded from this outcome denominator. |
| `vhd_valve_count_category` | categorical | 1  /  >=2  /  NA | Number of affected valves among detected VHD cases, classified as one or at least two from the supplied valve severity fields. Missing for participants without the primary event. |
| `vhd_aortic_stenosis_severity` | categorical | Mininal/Mild  /  Moderate  /  Severe  /  NA | Severity of aortic stenosis among detected VHD cases, based on the corresponding 90-day moderate and severe event records. Missing without the primary event. |
| `vhd_aortic_regurgitation_severity` | categorical | Mininal/Mild  /  Moderate  /  Severe  /  NA | Severity of aortic regurgitation among detected VHD cases, based on the corresponding 90-day moderate and severe event records. Missing without the primary event. |
| `vhd_mitral_regurgitation_severity` | categorical | Mininal/Mild  /  Moderate  /  Severe  /  NA | Severity of mitral regurgitation among detected VHD cases, based on the corresponding 90-day moderate and severe event records. Missing without the primary event. |
| `vhd_tricuspid_regurgitation_severity` | categorical | Mininal/Mild  /  Moderate  /  Severe  /  NA | Severity of tricuspid regurgitation among detected VHD cases, based on the corresponding 90-day moderate and severe event records. Missing without the primary event. |
| `vhd_lvef_percent` | numeric | Numeric %; blank = not applicable | LVEF among detected VHD cases. Missing without the primary event. |
| `vhd_lvef_category` | categorical | <=40%  /  41-50%  /  >50%  /  NA | Category based on vhd_lvef_percent: at most 40%, 41-50%, or above 50%. Missing without the primary event. |
| `vhd_lvedd_mm` | numeric | Numeric mm; blank = not applicable | Left ventricular end-diastolic diameter (LVEDD) among detected VHD cases. Missing without the primary event. |
| `vhd_lvedd_category` | categorical | <60 mm  /  >=60 mm  /  NA | Category based on vhd_lvedd_mm: below 60 mm or at least 60 mm. Missing without the primary event. |
| `vhd_pasp_mmhg` | numeric | Numeric mmHg; blank = not applicable | PASP among detected VHD cases. Missing without the primary event. |
| `vhd_pasp_category` | categorical | <50 mmHg  /  >=50 mmHg  /  NA | Category based on vhd_pasp_mmhg: below 50 mmHg or at least 50 mmHg. Missing without the primary event. |
| `vhd_postindex_dyspnea` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Documented dyspnea after VHD diagnosis and within 90 days after index time. Missing without the primary event. |
| `vhd_postindex_palpitations` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Documented palpitations after VHD diagnosis and within 90 days after index time. Missing without the primary event. |
| `vhd_postindex_angina_chest_pain` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Documented angina or chest pain after VHD diagnosis and within 90 days after index time. Missing without the primary event. |
| `vhd_postindex_syncope` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Documented syncope after VHD diagnosis and within 90 days after index time. Missing without the primary event. |
| `vhd_postindex_orthopnea` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Documented orthopnea after VHD diagnosis and within 90 days after index time. Missing without the primary event. |
| `vhd_postindex_edema` | binary | 0 = no qualifying condition or finding documented; 1 = qualifying condition or finding documented | Documented edema after VHD diagnosis and within 90 days after index time. Missing without the primary event. |
| `outcome_tte` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | TTE completed for any reason within the 90-day outcome window. |
| `tte_trigger` | categorical | Trial-triggered  /  Physician-triggered  /  NA | Supplied trigger of the follow-up TTE. Missing when no TTE is recorded. |
| `outcome_vhd_with_symptoms` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected moderate/severe VHD with at least one qualifying symptom documented after VHD diagnosis and within 90 days after index time, based on the six supplied post-index symptom flags. |
| `outcome_vhd_without_recorded_symptoms` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected moderate/severe VHD without a qualifying symptom documented after VHD diagnosis and within 90 days after index time. |
| `outcome_vhd_ntprobnp_ge900` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected moderate/severe VHD with NT-proBNP at least 900 pg/mL within 90 days after the index time. |
| `outcome_left_sided_vhd` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Moderate/severe aortic stenosis, aortic regurgitation or mitral regurgitation within 90 days. |
| `outcome_left_sided_vhd_lvedd_ge60` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected left-sided moderate/severe VHD with LVEDD at least 60 mm. |
| `outcome_severe_vhd` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Severe aortic stenosis, aortic regurgitation, mitral regurgitation or tricuspid regurgitation detected within 90 days. |
| `outcome_severe_vhd_surgical_criteria` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Severe VHD satisfying the study-defined surgical criteria: AS with LVEF at most 50% or NT-proBNP at least 900; AR with LVEF at most 50%, aortic root at least 55 mm, or bicuspid valve with root at least 50 mm; MR with AF, PASP at least 50 mmHg, left atrial diameter at least 55 mm, or moderate/severe TR; TR with TAPSE below 17 mm. For the MR criterion, AF may be documented at baseline or within 90 days after index time. |
| `outcome_vhd_any_examination` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with at least one of the five examination components, performed after VHD diagnosis and within 90 days after index time. |
| `outcome_vhd_chest_xray` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with a chest radiograph after VHD diagnosis and within 90 days after index time. |
| `outcome_vhd_ntprobnp_test` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with NT-proBNP testing after VHD diagnosis and within 90 days after index time. |
| `outcome_vhd_troponin_test` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with high-sensitivity troponin I testing after VHD diagnosis and within 90 days after index time. |
| `outcome_vhd_exercise_test` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with exercise testing after VHD diagnosis and within 90 days after index time. |
| `outcome_vhd_transesophageal_echo` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with transesophageal echocardiography after VHD diagnosis and within 90 days after index time. |
| `outcome_vhd_any_management` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with at least one management component, occurring after VHD diagnosis and within 90 days after index time. |
| `outcome_vhd_cardiology_referral` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with cardiology referral after VHD diagnosis and within 90 days after index time. Participants with regular baseline cardiology follow-up have missing/not-applicable values and are excluded only from this component denominator. |
| `outcome_vhd_any_medication_change` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with initiation or adjustment of heart-failure or antihypertensive medication after VHD diagnosis and within 90 days after index time. |
| `outcome_vhd_hf_medication_change` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with initiation or adjustment of ACEi/ARB, ARNI, beta blocker, MRA, SGLT2 inhibitor or loop diuretic after VHD diagnosis and within 90 days after index time. |
| `outcome_vhd_antihypertensive_change` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with initiation or adjustment of ACEi/ARB, beta blocker, calcium channel blocker, thiazide, alpha blocker or another antihypertensive after VHD diagnosis and within 90 days after index time. |
| `outcome_vhd_surgery_referral` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with cardiovascular surgery referral after VHD diagnosis and within 90 days after index time. Participants with regular baseline cardiovascular surgery follow-up have missing/not-applicable values and are excluded only from this component denominator. |
| `outcome_vhd_valve_intervention` | outcome | 0 = no qualifying event documented; 1 = qualifying event documented | Detected VHD with surgical or transcatheter valve intervention after VHD diagnosis and within 90 days after index time. |

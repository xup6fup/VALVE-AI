# Analysis input

The analysis reads the prepared, readable-name dataset directly from `data/VALVE_AI.csv`. An alternative CSV path may be supplied to `run_all.R` or through the `VALVE_AI_DATA` environment variable.

`VALVE_AI_template.csv` contains only the 114 column names. It contains no participant records and cannot reproduce the analyses. Definitions and coding are provided in the main README.

Participant-level data are not distributed in the public repository. The de-identified analysis dataset may be requested from the corresponding author, subject to approval following a simplified ethics review. Full EHR-linkable records and raw ECG signals require formal institutional ethics review and are not needed to run this project.

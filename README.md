# CODECHECK bundle — Contrast exposure and short-term AKI

This bundle reproduces the propensity-score overlap-weighting sensitivity analysis for:

**Acil Servis Hastalarında İntravenöz İyotlu Kontrast Maruziyeti ile Kısa Dönem Akut Böbrek Hasarı Arasındaki İlişki: Retrospektif Kohort Çalışması**

## Analysis
- Contrast-enhanced CT: n=244
- Noncontrast CT: n=124
- Total: n=368
- AKI: creatinine increase >=0.3 mg/dL or >=1.5x baseline in the 48–72 h window
- Propensity model: age, sex, baseline eGFR, urea, bicarbonate, grouped ED presentation, grouped CT region
- Weighting: overlap weights
- Outcome: overlap-weighted logistic regression
- SE: HC0-style robust/sandwich covariance for the weighted score equations, treating overlap weights as fixed in the outcome model

## Expected main result
OR 0.598; robust 95% CI 0.292–1.222; p=0.158.
Weighted AKI: contrast 9.67%; noncontrast 15.19%.

## Run
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python analysis.py
```

Outputs:
- outputs/overlap_weighting_results.csv
- outputs/covariate_balance.csv
- outputs/analysis_summary.txt

## Optional R cross-check
`optional_crosscheck/analysis_base_R.R` independently implements the same algorithm in base R.
R was not installed in the environment used to prepare this bundle, so that optional script has not yet been executed here.

## Privacy
`data/analysis_data.csv` is a minimized analysis extract. Original source IDs were removed and replaced by sequential study IDs; no dates or direct identifiers are included. This does not by itself authorize public release. Confirm ethics/institutional data-governance permission before public GitHub upload. If public release is not permitted, provide the dataset to the codechecker through an approved protected channel.


## Minimal public analysis dataset

This repository uses `data/analysis_minimal_full_reproducibility.csv`.

It contains only the variables required to refit the propensity-score model and reproduce the final overlap-weighted analysis:
- contrast exposure
- AKI outcome
- age
- female indicator
- baseline eGFR
- urea
- venous bicarbonate
- four-level presentation group
- four-level CT-region group

Original source record numbers, original row order, raw presentation codes, raw CT codes, creatinine values, other laboratory variables, dates, names, addresses, and direct identifiers are not included.

#!/usr/bin/env python3
from pathlib import Path
import csv, math
import numpy as np
import statsmodels.api as sm
from scipy.stats import norm

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "analysis_minimal_full_reproducibility.csv"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

with DATA.open(encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

X, t, y = [], [], []

for r in rows:
    cg = r["presentation_group"]
    ctg = r["ct_group"]
    X.append([
        1.0,
        float(r["age"]) / 10.0,
        int(r["female"]),
        float(r["egfr"]) / 10.0,
        float(r["urea"]) / 10.0,
        float(r["hco3"]) / 10.0,
        int(cg == "Abdominopelvic_GI_GU"),
        int(cg == "Trauma"),
        int(cg == "Other_systemic_neuro"),
        int(ctg == "Abd_pelvis"),
        int(ctg == "Head_neck"),
        int(ctg == "Multi_other"),
    ])
    t.append(int(r["exposure"]))
    y.append(int(r["aki"]))

X = np.asarray(X, dtype=float)
t = np.asarray(t, dtype=float)
y = np.asarray(y, dtype=float)

ps_model = sm.Logit(t, X).fit(disp=False)
ps = ps_model.predict(X)
ow = np.where(t == 1, 1 - ps, ps)

Z = np.column_stack([np.ones(len(t)), t])
fit = sm.GLM(y, Z, family=sm.families.Binomial(), freq_weights=ow).fit()

beta = np.asarray(fit.params, dtype=float)
mu = 1 / (1 + np.exp(-(Z @ beta)))
A = mu * (1 - mu)

bread = np.linalg.inv(Z.T @ ((ow * A)[:, None] * Z))
resid = y - mu
meat = Z.T @ (((ow ** 2) * (resid ** 2))[:, None] * Z)
V = bread @ meat @ bread
se = np.sqrt(np.diag(V))

b = beta[1]
se_b = se[1]
OR = math.exp(b)
lo = math.exp(b - 1.96 * se_b)
hi = math.exp(b + 1.96 * se_b)
p = 2 * (1 - norm.cdf(abs(b / se_b)))

def weighted_rate(group):
    mask = t == group
    return float(np.average(y[mask], weights=ow[mask]))

r_contrast = weighted_rate(1)
r_noncontrast = weighted_rate(0)

names = [
    "age10", "female", "egfr10", "urea10", "hco3_10",
    "compl_abd", "compl_trauma", "compl_other",
    "ct_abd", "ct_head", "ct_multi"
]
C = X[:, 1:]

def wmean(x, w):
    return np.sum(w * x) / np.sum(w)

def smd(x):
    xt = x[t == 1]
    xc = x[t == 0]
    den = math.sqrt((np.var(xt, ddof=1) + np.var(xc, ddof=1)) / 2)
    pre = (np.mean(xt) - np.mean(xc)) / den if den else 0.0
    post = (wmean(xt, ow[t == 1]) - wmean(xc, ow[t == 0])) / den if den else 0.0
    return pre, post

with (OUT / "overlap_weighting_results.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow([
        "analysis", "n", "or", "robust_se_log_or",
        "ci95_lower", "ci95_upper", "p_value",
        "weighted_aki_contrast", "weighted_aki_noncontrast"
    ])
    w.writerow([
        "propensity_score_overlap_weighting", len(t),
        OR, se_b, lo, hi, p, r_contrast, r_noncontrast
    ])

with (OUT / "covariate_balance.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["covariate", "smd_before", "smd_after"])
    for j, name in enumerate(names):
        w.writerow([name, *smd(C[:, j])])

(OUT / "analysis_summary.txt").write_text(
    f"N = {len(t)}\n"
    f"Contrast-enhanced CT = {int(np.sum(t == 1))}\n"
    f"Noncontrast CT = {int(np.sum(t == 0))}\n"
    f"OR = {OR:.6f}\n"
    f"Robust/sandwich SE(log OR) = {se_b:.6f}\n"
    f"95% CI = {lo:.6f} to {hi:.6f}\n"
    f"p = {p:.6f}\n"
    f"Weighted AKI contrast = {r_contrast:.6f}\n"
    f"Weighted AKI noncontrast = {r_noncontrast:.6f}\n",
    encoding="utf-8"
)

print(f"OR {OR:.3f}; robust 95% CI {lo:.3f}-{hi:.3f}; p={p:.3f}")

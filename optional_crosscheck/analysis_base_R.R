#!/usr/bin/env Rscript

d <- read.csv(
  "data/analysis_minimal_full_reproducibility.csv",
  stringsAsFactors = FALSE
)

X <- cbind(
  1,
  age10 = d$age / 10,
  female = d$female,
  egfr10 = d$egfr / 10,
  urea10 = d$urea / 10,
  hco3_10 = d$hco3 / 10,
  compl_abd = as.integer(d$presentation_group == "Abdominopelvic_GI_GU"),
  compl_trauma = as.integer(d$presentation_group == "Trauma"),
  compl_other = as.integer(d$presentation_group == "Other_systemic_neuro"),
  ct_abd = as.integer(d$ct_group == "Abd_pelvis"),
  ct_head = as.integer(d$ct_group == "Head_neck"),
  ct_multi = as.integer(d$ct_group == "Multi_other")
)

# Propensity-score model
psfit <- glm.fit(
  x = X,
  y = d$exposure,
  family = binomial()
)
ps <- psfit$fitted.values

# Overlap weights
ow <- ifelse(d$exposure == 1, 1 - ps, ps)

# Weighted logistic outcome model
Z <- cbind(1, d$exposure)
fit <- glm.fit(
  x = Z,
  y = d$aki,
  weights = ow,
  family = binomial()
)

beta <- fit$coefficients
mu <- plogis(as.vector(Z %*% beta))
A <- mu * (1 - mu)

# HC0-style robust/sandwich covariance
bread <- solve(t(Z) %*% (Z * as.vector(ow * A)))
resid <- d$aki - mu
meat <- t(Z) %*% (Z * as.vector((ow^2) * (resid^2)))
V <- bread %*% meat %*% bread
se <- sqrt(diag(V))

b <- beta[2]
seb <- se[2]

OR <- exp(b)
lo <- exp(b - 1.96 * seb)
hi <- exp(b + 1.96 * seb)
p <- 2 * pnorm(-abs(b / seb))

wrate <- function(g) {
  weighted.mean(
    d$aki[d$exposure == g],
    ow[d$exposure == g]
  )
}

cat(
  sprintf(
    "OR %.3f; robust 95%% CI %.3f-%.3f; p=%.3f\n",
    OR, lo, hi, p
  )
)

cat(
  sprintf(
    "Weighted AKI contrast %.4f; noncontrast %.4f\n",
    wrate(1), wrate(0)
  )
)

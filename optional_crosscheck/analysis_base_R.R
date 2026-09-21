#!/usr/bin/env Rscript
d <- read.csv("data/analysis_data.csv", stringsAsFactors=FALSE)
cardiopulmonary <- c(2,4,5,6,7,13,24)
abdominopelvic <- c(10,14,17,21,26,27)
trauma <- c(1,15,16,18,25)
complaint_group <- function(x) {
  if (x %in% cardiopulmonary) return("Cardiopulmonary")
  if (x %in% abdominopelvic) return("Abdominopelvic_GI_GU")
  if (x %in% trauma) return("Trauma")
  "Other_systemic_neuro"
}
ct_group <- function(ct1,ct2,ct3) {
  codes <- c(ct1,ct2,ct3); codes <- codes[codes != 0]; anatomical <- codes[codes != 10]
  if (length(anatomical)==0) return("Multi_other")
  g <- sapply(anatomical,function(x) {
    if (x %in% c(1,2,9)) return("Head_neck")
    if (x==3) return("Thorax")
    if (x %in% c(4,5)) return("Abd_pelvis")
    "Multi_other"
  })
  if ("Multi_other" %in% g || length(unique(g))>1) return("Multi_other")
  g[1]
}
d$complaint_group <- sapply(d$complaint1,complaint_group)
d$ct_group <- mapply(ct_group,d$ct1,d$ct2,d$ct3)
X <- cbind(1,age10=d$age/10,female=d$female,egfr10=d$egfr/10,urea10=d$urea/10,hco3_10=d$hco3/10,
           compl_abd=as.integer(d$complaint_group=="Abdominopelvic_GI_GU"),
           compl_trauma=as.integer(d$complaint_group=="Trauma"),
           compl_other=as.integer(d$complaint_group=="Other_systemic_neuro"),
           ct_abd=as.integer(d$ct_group=="Abd_pelvis"),ct_head=as.integer(d$ct_group=="Head_neck"),
           ct_multi=as.integer(d$ct_group=="Multi_other"))
psfit <- glm.fit(x=X,y=d$exposure,family=binomial())
ps <- psfit$fitted.values
ow <- ifelse(d$exposure==1,1-ps,ps)
Z <- cbind(1,d$exposure)
fit <- glm.fit(x=Z,y=d$aki,weights=ow,family=binomial())
beta <- fit$coefficients; mu <- plogis(as.vector(Z %*% beta)); A <- mu*(1-mu)
bread <- solve(t(Z) %*% (Z * as.vector(ow*A)))
resid <- d$aki-mu
meat <- t(Z) %*% (Z * as.vector((ow^2)*(resid^2)))
V <- bread %*% meat %*% bread
se <- sqrt(diag(V))
b <- beta[2]; seb <- se[2]
OR <- exp(b); lo <- exp(b-1.96*seb); hi <- exp(b+1.96*seb); p <- 2*pnorm(-abs(b/seb))
wr <- function(g) weighted.mean(d$aki[d$exposure==g],ow[d$exposure==g])
cat(sprintf("OR %.3f; robust 95%% CI %.3f-%.3f; p=%.3f\n",OR,lo,hi,p))
cat(sprintf("Weighted AKI contrast %.4f; noncontrast %.4f\n",wr(1),wr(0)))

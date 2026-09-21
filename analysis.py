#!/usr/bin/env python3
from pathlib import Path
import csv, math
import numpy as np
import statsmodels.api as sm
from scipy.stats import norm

ROOT=Path(__file__).resolve().parent
DATA=ROOT/"data"/"analysis_data.csv"
OUT=ROOT/"outputs"; OUT.mkdir(exist_ok=True)
with DATA.open(encoding="utf-8") as f:
    rows=list(csv.DictReader(f))

def fi(r,k): return int(float(r[k]))
def ff(r,k): return float(r[k])

cardiopulmonary={2,4,5,6,7,13,24}
abdominopelvic={10,14,17,21,26,27}
trauma={1,15,16,18,25}

def complaint_group(x):
    if x in cardiopulmonary: return "Cardiopulmonary"
    if x in abdominopelvic: return "Abdominopelvic_GI_GU"
    if x in trauma: return "Trauma"
    return "Other_systemic_neuro"

def ct_group(r):
    codes=[fi(r,"ct1"),fi(r,"ct2"),fi(r,"ct3")]
    codes=[x for x in codes if x!=0]
    anatomical=[x for x in codes if x!=10]
    if not anatomical: return "Multi_other"
    g=[]
    for x in anatomical:
        if x in {1,2,9}: g.append("Head_neck")
        elif x==3: g.append("Thorax")
        elif x in {4,5}: g.append("Abd_pelvis")
        else: g.append("Multi_other")
    if "Multi_other" in g or len(set(g))>1: return "Multi_other"
    return g[0]

X=[]; t=[]; y=[]
for r in rows:
    cg=complaint_group(fi(r,"complaint1")); ctg=ct_group(r)
    X.append([1,ff(r,"age")/10,fi(r,"female"),ff(r,"egfr")/10,ff(r,"urea")/10,ff(r,"hco3")/10,
              int(cg=="Abdominopelvic_GI_GU"),int(cg=="Trauma"),int(cg=="Other_systemic_neuro"),
              int(ctg=="Abd_pelvis"),int(ctg=="Head_neck"),int(ctg=="Multi_other")])
    t.append(fi(r,"exposure")); y.append(fi(r,"aki"))

X=np.asarray(X,float); t=np.asarray(t,float); y=np.asarray(y,float)
psfit=sm.Logit(t,X).fit(disp=False)
ps=psfit.predict(X)
ow=np.where(t==1,1-ps,ps)

Z=np.column_stack([np.ones(len(t)),t])
fit=sm.GLM(y,Z,family=sm.families.Binomial(),freq_weights=ow).fit()
beta=np.asarray(fit.params,float)
mu=1/(1+np.exp(-(Z@beta)))
A=mu*(1-mu)
bread=np.linalg.inv(Z.T@((ow*A)[:,None]*Z))
resid=y-mu
meat=Z.T@(((ow**2)*(resid**2))[:,None]*Z)
V=bread@meat@bread
se=np.sqrt(np.diag(V))

b=beta[1]; seb=se[1]
OR=math.exp(b); lo=math.exp(b-1.96*seb); hi=math.exp(b+1.96*seb)
p=2*(1-norm.cdf(abs(b/seb)))
wr=lambda g: float(np.average(y[t==g],weights=ow[t==g]))
r1=wr(1); r0=wr(0)

with (OUT/"overlap_weighting_results.csv").open("w",newline="",encoding="utf-8") as f:
    w=csv.writer(f); w.writerow(["analysis","n","or","robust_se_log_or","ci95_lower","ci95_upper","p_value","weighted_aki_contrast","weighted_aki_noncontrast"])
    w.writerow(["propensity_score_overlap_weighting",len(t),OR,seb,lo,hi,p,r1,r0])

names=["age10","female","egfr10","urea10","hco3_10","compl_abd","compl_trauma","compl_other","ct_abd","ct_head","ct_multi"]
C=X[:,1:]
def wmean(x,w): return np.sum(w*x)/np.sum(w)
def smd(x):
    xt=x[t==1]; xc=x[t==0]
    den=math.sqrt((np.var(xt,ddof=1)+np.var(xc,ddof=1))/2)
    return ((np.mean(xt)-np.mean(xc))/den,(wmean(xt,ow[t==1])-wmean(xc,ow[t==0]))/den)
with (OUT/"covariate_balance.csv").open("w",newline="",encoding="utf-8") as f:
    w=csv.writer(f); w.writerow(["covariate","smd_before","smd_after"])
    for j,n in enumerate(names): w.writerow([n,*smd(C[:,j])])

(OUT/"analysis_summary.txt").write_text(
    f"N = {len(t)}\nContrast-enhanced CT = {int(sum(t==1))}\nNoncontrast CT = {int(sum(t==0))}\n"
    f"OR = {OR:.6f}\nRobust/sandwich SE(log OR) = {seb:.6f}\n95% CI = {lo:.6f} to {hi:.6f}\n"
    f"p = {p:.6f}\nWeighted AKI contrast = {r1:.6f}\nWeighted AKI noncontrast = {r0:.6f}\n",
    encoding="utf-8")
print(f"OR {OR:.3f}; robust 95% CI {lo:.3f}-{hi:.3f}; p={p:.3f}")

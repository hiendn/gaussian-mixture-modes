#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=FALSE); f<-sub("^--file=","",args[grepl("^--file=",args)])
if(length(f)) setwd(dirname(dirname(normalizePath(f))))
if(dir.exists(".R-library")) .libPaths(c(normalizePath(".R-library"),.libPaths()))
source("R/core.R"); source("R/geometry.R"); source("R/exact_geometry.R")
g<-run_geometry(); exact<-run_exact_geometry()
stopifnot(exact$checks_passed==31)
expected_r<-c(.38225103987206115,.38225103987206115,.2680459329951420)
expected_v<-c(.29382574551738358,.29432001280708256,.71255941372803033)
stopifnot(max(abs(g$local_radius$table$local_radius-expected_r))<1e-13,
          max(abs(g$local_radius$table$local_variance-expected_v))<1e-13,
          abs(g$local_radius$table$common_variance[3]-1.0666729559241791)<1e-13)
fixed<-readRDS("results/fixed_budget.rds"); joint<-readRDS("results/joint_budget.rds")
X<-as.matrix(read.csv("data/fixed_budget.csv")); J<-as.matrix(read.csv("data/joint_budget.csv"))
stopifnot(nrow(X)==96,nrow(J)==12,ncol(X)==2,ncol(J)==2)
# Gaussian likelihoods in physical coordinates provide a different computation
# from the whitened complete-data implementation in R/fit.R.
physical_ll<-function(z,data) {
  Bi<-solve(z$B); covariance<-Bi%*%t(Bi); means<-t(solve(z$B,t(z$v)))
  mean(mixture_log_density(data,means,covariance,softmax(z$alpha)))
}
checks<-lapply(seq_along(c(fixed$steps,joint$steps)),function(i) {
  rec<-c(fixed$steps,joint$steps)[[i]]; dat<-if(i<=3)X else J
  old<-physical_ll(rec$anchor,dat); new<-physical_ll(rec$accepted,dat)
  penalty<-(if(i<=3).15 else .05)/2*sum(vapply(names(rec$anchor),function(nm)
    sum((rec$accepted[[nm]]-rec$anchor[[nm]])^2),0))
  stopifnot(abs(new-rec$log_likelihood_per_observation)<1e-12,
    abs(new-old-rec$log_likelihood_increase_per_observation)<1e-12,
    abs(penalty-rec$proximal_penalty)<1e-12,
    rec$exact_checks$minimum>=-1e-10,rec$inner_checks$maximum_violation<=1e-10,
    rec$log_likelihood_increase_per_observation>=rec$EM_auxiliary_increase_per_observation-1e-9,
    rec$EM_auxiliary_increase_per_observation>=penalty-1e-9)
  list(check=if(i<=3)paste0("fixed_",i) else "joint_1",physical_likelihood_error=abs(new-rec$log_likelihood_per_observation),
    penalty_error=abs(penalty-rec$proximal_penalty))
})
legacy<-jsonlite::fromJSON("legacy/python/operator/full_inner_results.json")
legacy_joint<-jsonlite::fromJSON("legacy/python/joint_budget/joint_budget_step_check.json")
gain<-vapply(fixed$steps,`[[`,0,"log_likelihood_increase_per_observation")
diff_fixed<-max(abs(gain-legacy$steps$log_likelihood_increase_per_observation))
diff_joint<-joint$steps[[1]]$log_likelihood_increase_per_observation-legacy_joint$likelihood_increase_per_observation
stopifnot(diff_fixed<5e-6,abs(diff_joint)<5e-5)
# Check the generated manuscript table against the saved numerical outputs.
tab<-read.csv("results/fitting_summary.csv")
manuscript_dir<-if(file.exists("manuscript/manuscript.tex"))"manuscript" else "."
latex<-paste(readLines(file.path(manuscript_dir,"generated/fitting_table.tex")),collapse="\n")
for(i in 1:4) for(nm in c("likelihood_gain","auxiliary_gain","proximal_penalty"))
  stopifnot(grepl(sprintf("%.8f",tab[[nm]][i]),latex,fixed=TRUE))
report<-list(passed=TRUE,exact_geometric_sign_checks=31,accepted_fitting_states=4,
  physical_coordinate_checks=checks,maximum_fixed_gain_difference_from_original_python=diff_fixed,
  joint_gain_difference_from_original_python=diff_joint,
  joint_penalty_difference_from_original_python=joint$steps[[1]]$proximal_penalty-legacy_joint$proximal_penalty,
  scope="Saved states and reported computations only; no global optimum or asymptotic convergence certificate.")
jsonlite::write_json(report,"verification/verification_results.json",auto_unbox=TRUE,pretty=TRUE,digits=17)
cat("PASS: 31 exact geometry signs, four saved fit states, physical Gaussian likelihoods, proximal penalties and generated manuscript table.\n")

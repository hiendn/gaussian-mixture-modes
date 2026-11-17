#!/usr/bin/env Rscript
# Retained from the independent inline R validation executed for this project.
# The density, constraint, inner-support, and ascent formulas below are the same
# as that executed check. The wrapper adds dependency/run-shape checks and a
# separate output filename, preserving the existing narrative review reports.
# Run from the repository root: Rscript verification/independent_fit_check.R
# This script reads saved states and data. It does not run an optimizer or
# source any production R helper functions.

if (dir.exists('.R-library')) .libPaths(c(normalizePath('.R-library'), .libPaths()))
if (!requireNamespace('jsonlite', quietly = TRUE)) stop('Install jsonlite to write the independent check report.')
fixed <- readRDS('results/fixed_budget.rds')
joint <- readRDS('results/joint_budget.rds')
Xfixed <- as.matrix(read.csv('data/fixed_budget.csv'))
Xjoint <- as.matrix(read.csv('data/joint_budget.csv'))

# Independent implementation. No production R functions are sourced.
lse <- function(y) { m <- max(y); log(sum(exp(y-m)))+m }
eigmin <- function(A) min(eigen((A+t(A))/2, symmetric=TRUE, only.values=TRUE)$values)
physical <- function(z, X) {
  precision <- crossprod(z$B)
  covariance <- solve(precision)
  mus <- t(vapply(seq_len(nrow(z$v)), function(j) solve(z$B,z$v[j,]), numeric(ncol(X))))
  logw <- z$alpha-lse(z$alpha)
  lc <- vapply(seq_len(nrow(mus)), function(j) {
    e <- sweep(X,2,mus[j,],'-')
    logw[j]-ncol(X)*log(2*pi)/2-as.numeric(determinant(covariance,logarithm=TRUE)$modulus)/2-rowSums((e%*%precision)*e)/2
  }, numeric(nrow(X)))
  ll <- mean(apply(lc,1,lse))
  R <- exp(lc-apply(lc,1,lse))
  list(ll=ll, R=R, covariance=covariance, mus=mus, component_log_density=lc)
}

inner_independent <- function(z, a, beta, joint) {
  p <- ncol(z$v); K <- nrow(z$v); I <- diag(p)
  pair <- do.call(rbind,lapply(1:K,function(j)cbind(j,setdiff(1:K,j))))
  out <- c(); budget <- numeric(K)
  for (k in seq_len(nrow(pair))) {
    j <- pair[k,1]; h <- pair[k,2]
    dt <- a$v[h,]-a$v[j,]; d0 <- sqrt(sum(dt^2)); dd <- z$v[h,]-z$v[j,]
    rt <- a$r[j]; rad <- z$r[j]
    s <- z$alpha[h]-z$alpha[j]; st <- a$alpha[h]-a$alpha[j]
    lbound <- (d0-rt)*(sum(dt*dd)/d0-rad)-(d0-rt)^2/2-rad^2/2
    out[paste0('odds_',k)] <- z$a[k]-exp(s-lbound)
    if (!joint) {
      out[paste0('separation_',k)] <- sum(dt*dd)/d0-sqrt(3.01+rad^2+2/rad^2)
    } else {
      nt <- a$nu[j]
      out[paste0('separation_',k)] <- 2*sum(dt*dd)-d0^2-(1.01+z$nu[j]+rad^2+z$nu[j]^2/(2*nt*rt^2)+nt*rt^2/(2*rad^4))
      out[paste0('trust_delta_',k)] <- .12-sqrt(sum((dd-dt)^2))
      out[paste0('trust_odds_',k)] <- .1-abs(s-st)
      P <- I-a$M[j,,]; dmd <- sum(dt*drop(P%*%dt)); E <- exp(st-d0^2/2+sqrt(2)*d0)
      grad_delta <- E*(2*drop(P%*%dt)+dmd*(sqrt(2)/d0-1)*dt)
      lo <- d0-.12; hi <- d0+.12; peak <- min(hi,max(lo,sqrt(2)))
      em <- exp(st+.1-peak^2/2+sqrt(2)*peak)
      G <- sqrt(1+(hi+sqrt(2))^2); H <- 1+sqrt(2)/lo
      L <- max(1e-8,em*(hi^2*(G^2+H)+2*G*sqrt(4*hi^2+hi^4)+2+2*hi))
      changeM <- z$M[j,,]-a$M[j,,]
      term <- E*dmd+sum(grad_delta*(dd-dt))+E*dmd*(s-st)-E*sum(tcrossprod(dt)*changeM)+L/2*(sum((dd-dt)^2)+(s-st)^2+sum(changeM^2))
      budget[j] <- budget[j]+term
    }
  }
  for (j in 1:K) {
    U <- matrix(0,p,p)
    for (k in which(pair[,1]==j)) {
      DL <- z$lam[k,,]-a$lam[k,,]; Da <- z$a[k]-a$a[k]
      # Equivalent independent form of the matrix majorant:
      # U=a*Lambda+(sqrt(kappa)*Da*I-DLambda/sqrt(kappa))^2/4.
      D <- 5*Da*I-DL/5
      U <- U+z$a[k]*z$lam[k,,]+D%*%D/4
    }
    out[paste0('first_matrix_',j)] <- eigmin(beta*I-tcrossprod(z$c[j,])-U)
    if (joint) out[paste0('budget_',j)] <- z$nu[j]-sum(diag(z$M[j,,]))-budget[j]/beta
  }
  out
}

exact_independent <- function(z, beta, joint) {
  p <- ncol(z$v); K <- nrow(z$v); I <- diag(p)
  pair <- do.call(rbind,lapply(1:K,function(j)cbind(j,setdiff(1:K,j))))
  nu <- if (joint) z$nu else rep(2,K)
  out <- c(); actual <- numeric(nrow(pair))
  for (k in seq_len(nrow(pair))) {
    j <- pair[k,1]; h <- pair[k,2]; dd <- z$v[h,]-z$v[j,]; d2 <- sum(dd^2)
    actual[k] <- exp(z$alpha[h]-z$alpha[j]-d2/2+z$r[j]*sqrt(d2))
    out[paste0('separation_',k)] <- d2-(nu[j]+1.01+z$r[j]^2+nu[j]/z$r[j]^2)
    out[paste0('odds_',k)] <- z$a[k]-actual[k]
    out[paste0('lambda_',k)] <- eigmin(z$lam[k,,])
    out[paste0('second_matrix_',k)] <- eigmin(beta*I+z$lam[k,,]-tcrossprod(dd-z$c[j,]))
  }
  for (j in 1:K) {
    first_actual <- first_lifted <- beta*I-tcrossprod(z$c[j,])
    budget <- 0
    for (k in which(pair[,1]==j)) {
      first_actual <- first_actual-actual[k]*z$lam[k,,]
      first_lifted <- first_lifted-z$a[k]*z$lam[k,,]
      h <- pair[k,2]; dd <- z$v[h,]-z$v[j,]; d2 <- sum(dd^2)
      M <- if (joint) z$M[j,,] else I
      budget <- budget+exp(z$alpha[h]-z$alpha[j]-d2/2+sqrt(2*d2))*drop(t(dd)%*%(I-M)%*%dd)
    }
    out[paste0('first_actual_',j)] <- eigmin(first_actual)
    out[paste0('first_lifted_',j)] <- eigmin(first_lifted)
    if (joint) {
      out[paste0('M_psd_',j)] <- eigmin(M)
      out[paste0('I_minus_M_psd_',j)] <- eigmin(I-M)
      out[paste0('budget_',j)] <- nu[j]-sum(diag(M))-budget/beta
      out[paste0('nu_lower_',j)] <- nu[j]-1
    }
  }
  if (joint) {
    domain <- c(diag(z$B)-.85,1.15-diag(z$B),z$B[1,2]+.15,.15-z$B[1,2],z$v+3,23-z$v,z$alpha+.25,.25-z$alpha,z$r-.02,1-z$r,z$nu-1,1e4-z$nu,z$a)
  } else {
    domain <- c(diag(z$B)-.3,3-diag(z$B),z$B[1,2]+2,2-z$B[1,2],z$v+8,8-z$v,z$alpha+3,3-z$alpha,z$r-.01,1-z$r,z$a)
  }
  out['domain'] <- min(domain)
  attr(out,'equality_residual') <- max(abs(c(z$B[2,1],z$alpha[K])))
  out
}

reports <- list(); records <- c(fixed$steps,joint$steps)
stopifnot(length(fixed$steps)==3L, length(joint$steps)==1L)
for (i in seq_along(records)) {
  rec <- records[[i]]; isjoint <- i==4
  X <- if (isjoint) Xjoint else Xfixed
  beta <- if (isjoint) .95 else .99; gamma <- if (isjoint) .05 else .15
  a <- rec$anchor; z <- rec$accepted
  old <- physical(a,X); new <- physical(z,X)
  # Physical complete Gaussian log densities include constants which cancel.
  qgain <- sum(old$R*(new$component_log_density-old$component_log_density))/nrow(X)
  penalty <- gamma/2*sum(vapply(names(a),function(nm)sum((z[[nm]]-a[[nm]])^2),numeric(1)))
  lgain <- new$ll-old$ll
  exact <- exact_independent(z,beta,isjoint)
  inner <- inner_independent(z,a,beta,isjoint)
  blend_error <- max(vapply(names(a),function(nm)max(abs(z[[nm]]-a[[nm]]-rec$retained_step_fraction*(rec$raw_trial[[nm]]-a[[nm]]))),numeric(1)))
  chain_error <- if (i>1 && i<4) max(vapply(names(a),function(nm)max(abs(a[[nm]]-records[[i-1]]$accepted[[nm]])),numeric(1))) else 0
  stored_errors <- c(log_likelihood=abs(new$ll-rec$log_likelihood_per_observation),likelihood_gain=abs(lgain-rec$log_likelihood_increase_per_observation),auxiliary_gain=abs(qgain-rec$EM_auxiliary_increase_per_observation),penalty=abs(penalty-rec$proximal_penalty))
  reports[[i]] <- list(check=if (isjoint) 'joint' else 'fixed',step=rec$step,solver_status=rec$solver_status,
    physical_log_likelihood=new$ll,likelihood_gain=lgain,auxiliary_gain=qgain,proximal_penalty=penalty,
    minimum_exact_slack=min(exact),minimum_inner_slack=min(inner),
    equality_residual=attr(exact,'equality_residual'),maximum_stored_value_discrepancy=max(stored_errors),
    stored_value_discrepancies=as.list(stored_errors),blend_error=blend_error,chain_error=chain_error,
    exact_slacks=as.list(exact),inner_slacks=as.list(inner),full_witness_fields=names(z),
    passed=min(exact)>-1e-9 && min(inner)>-1e-9 && attr(exact,'equality_residual')<1e-10 &&
      lgain>0 && lgain>=qgain-1e-12 && qgain>=penalty-1e-12 && max(stored_errors)<1e-12 &&
      blend_error<1e-12 && chain_error<1e-12)
}
original <- jsonlite::read_json('legacy/python/operator/full_inner_results.json',simplifyVector=TRUE)
oldjoint <- jsonlite::read_json('legacy/python/joint_budget/joint_budget_step_check.json',simplifyVector=TRUE)
refg <- original$steps$log_likelihood_increase_per_observation
comparison <- list(
  fixed_likelihood_gain_differences=vapply(reports[1:3],`[[`,0,'likelihood_gain')-refg,
  joint_likelihood_gain_difference=reports[[4]]$likelihood_gain-oldjoint$likelihood_increase_per_observation,
  joint_proximal_penalty_difference=reports[[4]]$proximal_penalty-oldjoint$proximal_penalty,
  interpretation='R and Python are separate conic solves of the same stated subproblems. The joint solve is optimal_inaccurate in both, and its whole-lift proximal penalty differs; exact feasibility and ascent are checked independently rather than asserting identical iterates.')
output <- list(
  scope='Independent computation review of three fixed-budget and one joint-budget accepted iterates. No production fitting helper functions sourced; direct physical Gaussian densities, independent nonlinear and inner-support formulas.',
  original_fixed_data_identical=identical(unname(Xfixed),unname(as.matrix(read.csv('legacy/python/operator/full_inner_data.csv')))),
  data_dimensions=list(fixed=dim(Xfixed),joint=dim(Xjoint)),
  covariance_discrepancies=list(fixed=max(abs(physical(fixed$final_lift,Xfixed)$covariance-fixed$final_shared_covariance)),joint=max(abs(physical(joint$final_lift,Xjoint)$covariance-joint$final_shared_covariance))),
  accepted_iterates=reports,comparison_with_original_python=comparison,
  passed=all(vapply(reports,`[[`,FALSE,'passed')),
  limits=c('Floating-point checks, not interval certificates of fitted states.','No independently verified convex primal-dual gap is available from this port.','No global nonconvex optimum, stationary-limit or asymptotic-convergence claim is established.'))
stopifnot(output$passed, output$original_fixed_data_identical)
jsonlite::write_json(output,'verification/independent_fit_check.json',pretty=TRUE,auto_unbox=TRUE,digits=17)
print(vapply(reports,function(r)c(LL=r$physical_log_likelihood,likelihood_gain=r$likelihood_gain,
  Q_gain=r$auxiliary_gain,penalty=r$proximal_penalty,minimum_exact=r$minimum_exact_slack,
  minimum_inner=r$minimum_inner_slack,stored_error=r$maximum_stored_value_discrepancy),numeric(7)))
cat('PASSED:',output$passed,'\n')

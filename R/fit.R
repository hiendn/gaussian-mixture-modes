# Native R convex-inner EM/MM implementation for the manuscript examples.
# CVXR >= 1.8 and Clarabel are required. All stored gaps are solver diagnostics,
# not gaps to the global nonconvex likelihood optimum.

fit_pairs <- function(K) do.call(rbind, lapply(seq_len(K), function(j)
  cbind(j = j, h = setdiff(seq_len(K), j))))
fit_lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
fit_eigmin <- function(A) min(eigen((A + t(A))/2, symmetric = TRUE, only.values = TRUE)$values)
fit_norm <- function(x) sqrt(sum(x^2))
fit_distance2 <- function(z, anchor) sum(vapply(names(z), function(nm)
  sum((z[[nm]] - anchor[[nm]])^2), numeric(1)))
fit_blend <- function(anchor, raw, fraction) setNames(lapply(names(anchor), function(nm)
  anchor[[nm]] + fraction * (raw[[nm]] - anchor[[nm]])), names(anchor))
fit_logkern <- function(z, X) {
  BX <- X %*% t(z$B)
  sapply(seq_len(nrow(z$v)), function(j)
    z$alpha[j] - rowSums(sweep(BX, 2, z$v[j, ], "-")^2)/2)
}
fit_ll <- function(z, X) mean(apply(fit_logkern(z, X), 1, fit_lse)) -
  fit_lse(z$alpha) + sum(log(diag(z$B))) - ncol(X)*log(2*pi)/2
fit_responsibilities <- function(z, X) {
  L <- fit_logkern(z, X)
  exp(L - apply(L, 1, fit_lse))
}
fit_aux <- function(z, X, R) {
  BX <- X %*% t(z$B)
  errors <- sapply(seq_len(nrow(z$v)), function(j)
    rowSums(sweep(BX, 2, z$v[j, ], "-")^2))
  sum(log(diag(z$B))) - sum(R * errors)/(2*nrow(X)) +
    sum(colMeans(R)*z$alpha) - fit_lse(z$alpha)
}

fit_initial_fixed <- function() {
  K <- 3L; p <- 2L; pairs <- fit_pairs(K)
  ang <- (0:2)*2*pi/3
  v <- 2.8/sqrt(3)*cbind(cos(ang), sin(ang))
  z <- list(B = diag(p), v = v, alpha = rep(0, K), r = rep(.68, K),
    c = matrix(0, K, p), a = rep(exp(-2.016), nrow(pairs)),
    lam = array(0, c(nrow(pairs), p, p)))
  for (j in seq_len(K)) {
    e1 <- -v[j, ]/fit_norm(v[j, ]); e2 <- c(-e1[2], e1[1])
    frame <- cbind(e1, e2); z$c[j, ] <- .51*e1
    for (k in which(pairs[, 1] == j)) {
      sgn <- sign(sum((v[pairs[k, 2], ]-v[j, ])*e2))
      z$lam[k, , ] <- frame %*% matrix(c(2.72, sgn*2.70, sgn*2.70, 2.70), 2) %*% t(frame)
    }
  }
  z
}

fit_initial_joint <- function() {
  K <- 3L; p <- 2L; pairs <- fit_pairs(K)
  z <- list(B = diag(p), v = rbind(c(0,0), c(2.3,0), c(0,20)), alpha = rep(0,K),
    r = c(.6,.6,.06), nu = rep(1.001,K), M = array(0,c(K,p,p)),
    c = matrix(0,K,p), a = rep(0,nrow(pairs)), lam = array(0,c(nrow(pairs),p,p)))
  z$M[1,,] <- z$M[2,,] <- diag(c(1,0)); z$M[3,,] <- diag(c(0,1))
  for (j in seq_len(K)) {
    ids <- which(pairs[,1] == j)
    deltas <- sweep(z$v[pairs[ids,2],,drop=FALSE], 2, z$v[j,], "-")
    ds <- sqrt(rowSums(deltas^2)); caps <- exp(-ds^2/2+z$r[j]*ds)
    centre <- colSums(deltas*caps)/(1+sum(caps)); z$c[j,] <- centre
    for (i in seq_along(ids)) {
      k <- ids[i]; z$a[k] <- caps[i]
      z$lam[k,,] <- (sum((deltas[i,]-centre)^2)-.95+1e-4)*diag(p)
    }
  }
  z
}
fit_joint_data <- function() {
  offsets <- rbind(c(-.20,.10),c(.15,-.20),c(.10,.25),c(-.05,-.10))
  centres <- rbind(c(0,0),c(2.3,0),c(0,20))+rbind(c(.02,0),c(.04,.05),c(.12,-.02))
  do.call(rbind,lapply(seq_len(3), function(j) sweep(offsets,2,centres[j,],"+")))
}

# Every returned constraint is a signed slack, feasible when nonnegative.
fit_exact_checks <- function(z, beta, tau = .01, fixed_nu = 2) {
  K <- nrow(z$v); p <- ncol(z$v); I <- diag(p); pairs <- fit_pairs(K)
  nu <- if (is.null(z$nu)) rep(fixed_nu,K) else z$nu
  names_pair <- c("separation","odds_epigraph","lambda_psd","second_matrix")
  pair <- matrix(NA_real_,nrow(pairs),length(names_pair),dimnames=list(NULL,names_pair))
  A <- numeric(nrow(pairs))
  for (k in seq_len(nrow(pairs))) {
    j <- pairs[k,1]; h <- pairs[k,2]; delta <- z$v[h,]-z$v[j,]; d <- fit_norm(delta)
    A[k] <- exp(z$alpha[h]-z$alpha[j]-d*d/2+z$r[j]*d)
    dc <- delta-z$c[j,]
    pair[k,] <- c(d*d-(1+tau+nu[j]+z$r[j]^2+nu[j]/z$r[j]^2),z$a[k]-A[k],
      fit_eigmin(z$lam[k,,]),fit_eigmin(beta*I+z$lam[k,,]-tcrossprod(dc)))
  }
  comp_names <- c("first_lifted_matrix","first_actual_matrix")
  if (!is.null(z$M)) comp_names <- c(comp_names,"M_psd","I_minus_M_psd","budget","nu_lower")
  component <- matrix(NA_real_,K,length(comp_names),dimnames=list(NULL,comp_names))
  for (j in seq_len(K)) {
    ids <- which(pairs[,1]==j)
    lift <- Reduce(`+`,lapply(ids,function(k) z$a[k]*z$lam[k,,]))
    actual <- Reduce(`+`,lapply(ids,function(k) A[k]*z$lam[k,,]))
    component[j,1:2] <- c(fit_eigmin(beta*I-tcrossprod(z$c[j,])-lift),
      fit_eigmin(beta*I-tcrossprod(z$c[j,])-actual))
    if (!is.null(z$M)) {
      residual <- sum(vapply(ids,function(k) {
        h <- pairs[k,2]; delta <- z$v[h,]-z$v[j,]; d <- fit_norm(delta)
        exp(z$alpha[h]-z$alpha[j]-d*d/2+sqrt(min(p,K-1))*d)*
          drop(t(delta) %*% (I-z$M[j,,]) %*% delta)
      },numeric(1)))
      component[j,3:6] <- c(fit_eigmin(z$M[j,,]),fit_eigmin(I-z$M[j,,]),
        nu[j]-sum(diag(z$M[j,,]))-residual/beta,nu[j]-1)
    }
  }
  minima <- c(apply(pair,2,min),apply(component,2,min))
  list(minima=as.list(minima), pairs=cbind(pairs,pair), components=component,
    minimum=min(minima), weights=exp(z$alpha-fit_lse(z$alpha)), radii=z$r)
}

fit_curvature <- function(d, s, rho=.12, rho_s=.1, q=2, floor=1e-8) {
  lower <- d-rho; upper <- d+rho
  if (lower <= 0) stop("The displacement trust ball must avoid zero.")
  peak <- min(max(sqrt(q),lower),upper)
  Amax <- exp(s+rho_s-peak^2/2+sqrt(q)*peak)
  G <- sqrt(1+(upper+sqrt(q))^2); H <- 1+sqrt(q)/lower
  max(floor,Amax*(upper^2*(G^2+H)+2*G*sqrt(4*upper^2+upper^4)+2+2*upper))
}

fit_inner_checks <- function(z, anchor, beta, joint=FALSE, tau=.01, kappa=25,
                             rho=.12, rho_s=.1, fixed_nu=2) {
  p <- ncol(z$v); K <- nrow(z$v); I <- diag(p); pairs <- fit_pairs(K)
  pair <- matrix(NA_real_,nrow(pairs),if(joint) 4 else 2)
  colnames(pair) <- if(joint) c("separation","odds_epigraph","displacement_trust","log_odds_trust") else c("separation","odds_epigraph")
  terms <- vector("list",K)
  for (k in seq_len(nrow(pairs))) {
    j <- pairs[k,1]; h <- pairs[k,2]; dt <- anchor$v[h,]-anchor$v[j,]; d0 <- fit_norm(dt)
    delta <- z$v[h,]-z$v[j,]; rad <- z$r[j]; tangent <- d0-anchor$r[j]
    low <- tangent*(sum(dt/d0*delta)-rad)-tangent^2/2-rad^2/2
    pair[k,2] <- z$a[k]-exp(z$alpha[h]-z$alpha[j]-low)
    if (!joint) pair[k,1] <- sum(dt/d0*delta)-sqrt(fixed_nu+1+tau+rad^2+fixed_nu/rad^2)
    else {
      nt <- anchor$nu[j]; rt <- anchor$r[j]; s <- z$alpha[h]-z$alpha[j]
      st <- anchor$alpha[h]-anchor$alpha[j]
      sep <- 1+tau+z$nu[j]+rad^2+z$nu[j]^2/(2*nt*rt^2)+nt*rt^2*rad^(-4)/2
      pair[k,1] <- 2*sum(dt*delta)-d0^2-sep
      pair[k,3:4] <- c(rho-fit_norm(delta-dt),rho_s-abs(s-st))
      Pm <- I-anchor$M[j,,]; b0 <- drop(t(dt)%*%Pm%*%dt); A0 <- exp(st-d0^2/2+sqrt(2)*d0)
      gd <- A0*(2*drop(Pm%*%dt)+b0*(-dt+sqrt(2)*dt/d0))
      dM <- z$M[j,,]-anchor$M[j,,]; L <- fit_curvature(d0,st,rho,rho_s)
      term <- A0*b0+sum(gd*(delta-dt))+A0*b0*(s-st)-A0*sum(tcrossprod(dt)*dM)+
        L*(sum((delta-dt)^2)+(s-st)^2+sum(dM^2))/2
      terms[[j]] <- c(terms[[j]],term)
    }
  }
  cm <- matrix(NA_real_,K,if(joint)2 else 1)
  colnames(cm) <- if(joint)c("first_matrix","budget") else "first_matrix"
  for (j in seq_len(K)) {
    U <- matrix(0,p,p)
    for(k in which(pairs[,1]==j)) {
      at <- anchor$a[k]; Lt <- anchor$lam[k,,]; aa <- z$a[k]; LL <- z$lam[k,,]
      W <- sqrt(kappa)*(aa-at)*I+(LL-Lt)/sqrt(kappa)
      U <- U+at*LL+aa*Lt-at*Lt+W%*%W/4
    }
    cm[j,1] <- fit_eigmin(beta*I-tcrossprod(z$c[j,])-U)
    if(joint) cm[j,2] <- z$nu[j]-sum(diag(z$M[j,,]))-sum(terms[[j]])/beta
  }
  # Domain restrictions from the corresponding exact convex subproblem.
  if(joint) domain <- c(diag(z$B)-.85,1.15-diag(z$B),z$B[1,2]+.15,.15-z$B[1,2],
    z$v+3,23-z$v,z$alpha+.25,.25-z$alpha,z$r-.02,1-z$r,z$nu-1,1e4-z$nu,z$a)
  else domain <- c(diag(z$B)-.3,3-diag(z$B),z$B[1,2]+2,2-z$B[1,2],
    z$v+8,8-z$v,z$alpha+3,3-z$alpha,z$r-.01,1-z$r,z$a)
  equality <- max(abs(c(z$B[2,1],z$alpha[K])))
  minima <- c(apply(pair,2,min),apply(cm,2,min),domain=min(domain),equalities=-equality)
  list(minima=as.list(minima),pairs=cbind(pairs,pair),components=cm,
    maximum_violation=max(0,-min(minima)))
}

fit_require_solver <- function() {
  if (!requireNamespace("CVXR",quietly=TRUE) || utils::packageVersion("CVXR") < "1.8")
    stop("Install CVXR >= 1.8 to run the native R conic fitting checks.")
  if (!requireNamespace("clarabel",quietly=TRUE)) stop("Install the R clarabel package.")
}
fit_sum <- function(xs) Reduce(`+`,xs)
fit_var <- function(n,m=1L,...) CVXR::Variable(c(n,m),...)
fit_dot <- function(a,b) sum(CVXR::multiply(matrix(a,ncol=1),b))
fit_psd <- function(A) CVXR::PSD(A)

fit_solve_step <- function(anchor,X,joint=FALSE,beta=if(joint).95 else .99,
                           gamma=if(joint).05 else .15,kappa=25,tau=.01,
                           rho=.12,rho_s=.1,tolerance=1e-9,max_iter=250,
                           acceptance_tolerance=if(joint)1e-7 else 1e-10) {
  fit_require_solver()
  p <- ncol(X); K <- nrow(anchor$v)
  if(p!=2L || K!=3L) stop("This port currently implements the manuscript's p=2, K=3 envelopes.")
  pairs <- fit_pairs(K); m <- nrow(pairs); I <- diag(p); n <- nrow(X)
  B <- fit_var(p,p,name="B"); v <- lapply(seq_len(K),function(j)fit_var(p,name=paste0("v",j)))
  al <- fit_var(K,name="alpha"); r <- fit_var(K,nonneg=TRUE,name="r")
  cc <- lapply(seq_len(K),function(j)fit_var(p,name=paste0("c",j)))
  a <- fit_var(m,nonneg=TRUE,name="a")
  lam <- lapply(seq_len(m),function(k)fit_var(p,p,symmetric=TRUE,name=paste0("Lambda",k)))
  if(joint) {
    nu <- fit_var(K,name="nu"); Ms <- lapply(seq_len(K),function(j)fit_var(p,p,symmetric=TRUE,name=paste0("M",j)))
    constraints <- list(B[2,1]==0,CVXR::DiagMat(B)>=.85,CVXR::DiagMat(B)<=1.15,B[1,2]>=-.15,B[1,2]<=.15,
      al>=-.25,al<=.25,al[K]==0,r>=.02,r<=1,nu>=1,nu<=1e4)
    for(j in seq_len(K)) constraints <- c(constraints,list(v[[j]]>=-3,v[[j]]<=23,fit_psd(Ms[[j]]),fit_psd(I-Ms[[j]])))
  } else {
    constraints <- list(B[2,1]==0,CVXR::DiagMat(B)>=.3,CVXR::DiagMat(B)<=3,B[1,2]>=-2,B[1,2]<=2,
      al>=-3,al<=3,al[K]==0,r>=.01,r<=1)
    for(j in seq_len(K)) constraints <- c(constraints,list(v[[j]]>=-8,v[[j]]<=8))
  }
  budget <- vector("list",K); curvatures <- numeric(m)
  for(k in seq_len(m)) {
    j <- pairs[k,1]; h <- pairs[k,2]; dt <- anchor$v[h,]-anchor$v[j,]; d0 <- fit_norm(dt)
    delta <- v[[h]]-v[[j]]; s <- al[h]-al[j]; s0 <- d0-anchor$r[j]
    low <- s0*(fit_dot(dt/d0,delta)-r[j])-s0^2/2-CVXR::square(r[j])/2
    constraints <- c(constraints,list(exp(s-low)<=a[k],fit_psd(lam[[k]])))
    if(joint) {
      rt <- anchor$r[j]; nt <- anchor$nu[j]; st <- anchor$alpha[h]-anchor$alpha[j]
      sep <- 1+tau+nu[j]+CVXR::square(r[j])+CVXR::square(nu[j])/(2*nt*rt^2)+
        nt*rt^2*CVXR::power(r[j],-4)/2
      constraints <- c(constraints,list(sep<=2*fit_dot(dt,delta)-d0^2,
        CVXR::p_norm(delta-matrix(dt,ncol=1),2)<=rho,abs(s-st)<=rho_s))
      Pm <- I-anchor$M[j,,]; b0 <- drop(t(dt)%*%Pm%*%dt); A0 <- exp(st-d0^2/2+sqrt(2)*d0)
      gd <- A0*(2*drop(Pm%*%dt)+b0*(-dt+sqrt(2)*dt/d0)); dM <- Ms[[j]]-anchor$M[j,,]
      L <- fit_curvature(d0,st,rho,rho_s); curvatures[k] <- L
      term <- A0*b0+fit_dot(gd,delta-matrix(dt,ncol=1))+A0*b0*(s-st)-
        A0*sum(CVXR::multiply(tcrossprod(dt),dM))+L/2*(
          CVXR::sum_squares(delta-matrix(dt,ncol=1))+CVXR::square(s-st)+CVXR::sum_squares(dM))
      budget[[j]] <- c(budget[[j]],list(term))
    } else {
      H <- CVXR::p_norm(CVXR::vstack(sqrt(3+tau),r[j],sqrt(2)*CVXR::inv_pos(r[j])),2)
      constraints <- c(constraints,list(H<=fit_dot(dt/d0,delta)))
    }
    dc <- delta-cc[[j]]
    constraints <- c(constraints,list(fit_psd(CVXR::bmat(list(
      list(beta*I+lam[[k]],dc),list(t(dc),matrix(1,1,1)))))))
  }
  for(j in seq_len(K)) {
    if(joint) constraints <- c(constraints,list(sum(CVXR::DiagMat(Ms[[j]]))+fit_sum(budget[[j]])/beta<=nu[j]))
    linear <- matrix(0,p,p); Ws <- list()
    for(k in which(pairs[,1]==j)) {
      at <- anchor$a[k]; Lt <- anchor$lam[k,,]
      linear <- linear+at*lam[[k]]+a[k]*Lt-at*Lt
      Ws <- c(Ws,list((sqrt(kappa)*(a[k]-at)*I+(lam[[k]]-Lt)/sqrt(kappa))/2))
    }
    blocks <- list(c(list(beta*I-linear,cc[[j]]),Ws),
      c(list(t(cc[[j]]),matrix(1,1,1)),lapply(Ws,function(x)matrix(0,1,p))))
    for(i in seq_along(Ws)) blocks <- c(blocks,list(c(list(Ws[[i]],matrix(0,p,1)),
      lapply(seq_along(Ws),function(h)if(i==h)I else matrix(0,p,p)))))
    constraints <- c(constraints,list(fit_psd(CVXR::bmat(blocks))))
  }
  R <- fit_responsibilities(anchor,X); BX <- X%*%t(B)
  error_terms <- lapply(seq_len(K),function(j) {
    difference <- BX-matrix(1,n,1)%*%t(v[[j]])
    CVXR::sum_squares(CVXR::multiply(matrix(sqrt(R[,j]),n,p),difference))
  })
  Qt <- sum(log(CVXR::DiagMat(B)))-fit_sum(error_terms)/(2*n)+fit_dot(colMeans(R),al)-CVXR::log_sum_exp(al)
  prox <- list(CVXR::sum_squares(B-anchor$B),CVXR::sum_squares(al-matrix(anchor$alpha,ncol=1)),
    CVXR::sum_squares(r-matrix(anchor$r,ncol=1)),CVXR::sum_squares(a-matrix(anchor$a,ncol=1)))
  for(j in seq_len(K)) prox <- c(prox,list(CVXR::sum_squares(v[[j]]-matrix(anchor$v[j,],ncol=1)),
    CVXR::sum_squares(cc[[j]]-matrix(anchor$c[j,],ncol=1))))
  for(k in seq_len(m)) prox <- c(prox,list(CVXR::sum_squares(lam[[k]]-anchor$lam[k,,])))
  if(joint) {
    prox <- c(prox,list(CVXR::sum_squares(nu-matrix(anchor$nu,ncol=1))))
    for(j in seq_len(K)) prox <- c(prox,list(CVXR::sum_squares(Ms[[j]]-anchor$M[j,,])))
  }
  problem <- CVXR::Problem(CVXR::Maximize(Qt-gamma*fit_sum(prox)/2),constraints)
  if(!CVXR::is_dcp(problem)) stop("The constructed R subproblem is not DCP.")
  started <- proc.time()[[3]]
  opt <- CVXR::psolve(problem,solver="CLARABEL",verbose=FALSE,tol_gap_abs=tolerance,
    tol_gap_rel=tolerance,tol_feas=tolerance,max_iter=max_iter)
  elapsed <- proc.time()[[3]]-started; solver_status <- CVXR::status(problem)
  if(!solver_status %in% c("optimal","optimal_inaccurate")) stop(paste("Solver status",solver_status))
  get_matrix <- function(x) as.matrix(CVXR::value(x)); get_vector <- function(x) as.numeric(CVXR::value(x))
  raw <- list(B=get_matrix(B),v=do.call(rbind,lapply(v,get_vector)),alpha=get_vector(al),r=get_vector(r),
    c=do.call(rbind,lapply(cc,get_vector)),a=get_vector(a),lam=array(0,c(m,p,p)))
  for(k in seq_len(m)) raw$lam[k,,] <- get_matrix(lam[[k]])
  if(joint) {
    raw$nu <- get_vector(nu); raw$M <- array(0,c(K,p,p))
    for(j in seq_len(K)) raw$M[j,,] <- get_matrix(Ms[[j]])
    # Preserve the same field order in the anchor and every stored state.
    raw <- raw[names(anchor)]
  }
  fraction <- .999
  for(attempt in seq_len(13)) {
    accepted <- fit_blend(anchor,raw,fraction)
    exact <- fit_exact_checks(accepted,beta,tau)
    inner <- fit_inner_checks(accepted,anchor,beta,joint,tau,kappa,rho,rho_s)
    qgain <- fit_aux(accepted,X,R)-fit_aux(anchor,X,R)
    penalty <- gamma*fit_distance2(accepted,anchor)/2; lgain <- fit_ll(accepted,X)-fit_ll(anchor,X)
    feasible <- exact$minimum>=-acceptance_tolerance && inner$maximum_violation<=acceptance_tolerance
    ascent <- qgain>=penalty-1e-9 && lgain>=qgain-1e-9 && lgain>0
    if(feasible && ascent) break
    fraction <- fraction*.8
  }
  if(!feasible || !ascent) stop("No retained trial passed the direct feasibility and ascent checks.")
  stats <- CVXR::solver_stats(problem)
  # Solver statistics APIs differ. Keep an auditable text representation and
  # expose no fabricated primal-dual or global-optimum gap.
  stats_text <- capture.output(str(stats))
  list(state=accepted,record=list(solver_status=solver_status,solver="CLARABEL",
    CVXR_version=as.character(utils::packageVersion("CVXR")),
    clarabel_version=as.character(utils::packageVersion("clarabel")),
    dcp_verified=TRUE,solve_wall_seconds=elapsed,subproblem_objective=as.numeric(opt),
    solver_statistics=stats_text,solver_gap=list(available=FALSE,
      reason="No independently verified primal-dual gap is extracted by this port; solver status and tolerances are retained."),
    tolerance=tolerance,max_iter=max_iter,retained_step_fraction=fraction,
    acceptance_attempts=attempt,log_likelihood_per_observation=fit_ll(accepted,X),
    log_likelihood_increase_per_observation=lgain,EM_auxiliary_increase_per_observation=qgain,
    proximal_penalty=penalty,exact_checks=exact,inner_checks=inner,
    pair_curvature_bounds=if(joint)curvatures else NULL,
    anchor=anchor,raw_trial=raw,accepted=accepted))
}

fit_fixed <- function(X,steps=3L,...) {
  X <- as.matrix(X); z <- fit_initial_fixed(); initial <- z; initial_ll <- fit_ll(z,X)
  records <- vector("list",steps)
  for(t in seq_len(steps)) {
    run <- fit_solve_step(z,X,joint=FALSE,...); z <- run$state
    records[[t]] <- c(list(step=t),run$record)
    message(sprintf("Fixed budget step %d: %s, likelihood gain %.10f",t,
      run$record$solver_status,run$record$log_likelihood_increase_per_observation))
  }
  Bi <- solve(z$B)
  list(scope="Native R fixed-budget bounded fitting check; no stationary or global optimum diagnosis.",
    settings=list(n=nrow(X),p=ncol(X),K=3,q=2,nu=2,beta=.99,tau=.01,
      proximal_gamma=.15,matrix_majorant_scale=25,fitting_steps=steps),
    initial_state=initial,initial_log_likelihood_per_observation=initial_ll,
    steps=records,total_log_likelihood_increase_per_observation=fit_ll(z,X)-initial_ll,
    final_shared_covariance=Bi%*%t(Bi),final_original_coordinate_means=t(solve(z$B,t(z$v))),final_lift=z)
}

fit_joint <- function(X=fit_joint_data(),...) {
  z <- fit_initial_joint(); run <- fit_solve_step(z,as.matrix(X),joint=TRUE,...)
  Bi <- solve(run$state$B)
  list(scope="Native R joint-budget one-step check; no asymptotic convergence or global optimum claim.",
    settings=list(n=nrow(X),p=ncol(X),K=3,q=2,beta=.95,tau=.01,proximal_gamma=.05,
      matrix_majorant_scale=25,displacement_trust_radius=.12,log_odds_trust_radius=.1),
    initial_exact_checks=fit_exact_checks(z,.95),steps=list(c(list(step=1),run$record)),
    final_lift=run$state,final_shared_covariance=Bi%*%t(Bi))
}

fit_fixed_density_sdp <- function() {
  fit_require_solver(); p <- 2L; I <- diag(p); d <- 2.8; radius <- .68
  cap <- exp(-d*d/2+d*radius); deltas <- rbind(c(sqrt(3)*d/2,d/2),c(sqrt(3)*d/2,-d/2))
  cc <- fit_var(p); level <- fit_var(1); L <- lapply(seq_len(2),function(h)fit_var(p,p,symmetric=TRUE))
  cons <- list(level>=0,fit_psd(CVXR::bmat(list(list(level*I-cap*(L[[1]]+L[[2]]),cc),
    list(t(cc),matrix(1,1,1))))))
  for(h in seq_len(2)) {
    dc <- matrix(deltas[h,],ncol=1)-cc
    cons <- c(cons,list(fit_psd(L[[h]]),fit_psd(CVXR::bmat(list(list(level*I+L[[h]],dc),
      list(t(dc),matrix(1,1,1)))))))
  }
  prob <- CVXR::Problem(CVXR::Minimize(level),cons)
  value <- CVXR::psolve(prob,solver="CLARABEL",tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10,max_iter=250)
  analytic <- (2*cap)*(sqrt(3)*d/2)^2/(1+2*cap)^2
  if(abs(value-analytic)>2e-8) stop("Fixed-density SDP disagrees with the analytic value.")
  centre <- as.numeric(CVXR::value(cc))
  multipliers <- lapply(L,function(x)as.matrix(CVXR::value(x)))
  first_slack <- fit_eigmin(as.numeric(value)*I-cap*(multipliers[[1]]+multipliers[[2]])-tcrossprod(centre))
  second_slacks <- vapply(seq_len(2),function(h)fit_eigmin(as.numeric(value)*I+multipliers[[h]]-tcrossprod(deltas[h,]-centre)),numeric(1))
  list(solver_status=CVXR::status(prob),sdp_value=as.numeric(value),analytic_value=analytic,
    absolute_difference=abs(value-analytic),radius=radius,centre=centre,multipliers=multipliers,
    first_matrix_slack=first_slack,second_matrix_slacks=second_slacks,
    checks_are_floating_point=TRUE)
}

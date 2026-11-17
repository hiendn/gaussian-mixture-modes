# Shared-covariance Gaussian calculations. All rows of means are components.
log_sum_exp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
softmax <- function(x) exp(x - log_sum_exp(x))
min_eigenvalue <- function(A) min(eigen((A + t(A))/2, symmetric=TRUE, only.values=TRUE)$values)
rho_squared <- function(d, nu, tau=0) {
  a <- d*d - nu - 1 - tau
  if (any(a < 2*sqrt(nu))) stop("Separation does not permit a real radius")
  2*nu/(a + sqrt(pmax(0, a*a - 4*nu)))
}
mixture_log_density <- function(X, means, covariance=diag(ncol(means)), weights=rep(1/nrow(means),nrow(means))) {
  X <- matrix(X,ncol=ncol(means)); U <- chol(covariance)
  vals <- vapply(seq_len(nrow(means)),function(j) {
    residual <- sweep(X,2,means[j,],"-")
    whitened <- forwardsolve(t(U),t(residual))
    log(weights[j]) - ncol(means)*log(2*pi)/2 - sum(log(diag(U))) - colSums(whitened^2)/2
  }, numeric(nrow(X)))
  apply(matrix(vals,nrow=nrow(X)),1,log_sum_exp)
}
posterior_summary <- function(y,means,weights=rep(1/nrow(means),nrow(means))) {
  residual <- sweep(means,2,y,"-")
  w <- softmax(log(weights)-rowSums(residual^2)/2)
  m <- colSums(means*w); centered <- sweep(means,2,m,"-")
  list(mean=m, covariance=crossprod(centered,centered*w),weights=w)
}
locate_modes <- function(means,weights=rep(1/nrow(means),nrow(means)),tol=1e-14,maxit=10000L) {
  ans <- lapply(seq_len(nrow(means)),function(j) {
    y <- means[j,]; converged <- FALSE
    for(i in seq_len(maxit)) {
      new <- posterior_summary(y,means,weights)$mean
      if(sqrt(sum((new-y)^2))<tol) { y<-new; converged<-TRUE; break }; y<-new
    }
    s <- posterior_summary(y,means,weights)
    if(!converged) stop("Posterior-mean iteration did not converge")
    list(mode=y,residual=sqrt(sum((s$mean-y)^2)),iterations=i,
         maximum_log_density_hessian_eigenvalue=max(eigen(s$covariance-diag(ncol(means)),symmetric=TRUE)$values))
  })
  list(modes=do.call(rbind,lapply(ans,`[[`,"mode")),residuals=vapply(ans,`[[`,0,"residual"),
       iterations=vapply(ans,`[[`,0L,"iterations"),
       maximum_log_density_hessian_eigenvalues=vapply(ans,`[[`,0,"maximum_log_density_hessian_eigenvalue"),
       scope="Modes located from component means only; exhaustiveness is supplied by the manuscript certificate.")
}
saturated_moments <- function(means,weights,j,radius) {
  means<-as.matrix(means); h<-setdiff(seq_len(nrow(means)),j)
  delta<-sweep(means[h,,drop=FALSE],2,means[j,],"-"); d2<-rowSums(delta^2)
  caps<-weights[h]/weights[j]*exp(-d2/2+radius*sqrt(d2)); S<-sum(caps)
  m<-colSums(delta*caps)/(1+S)
  C<-crossprod(delta,delta*caps)/(1+S)-tcrossprod(m)
  list(odds=caps,S=S,mean=m,covariance=C,trace=sum(diag(C)),eigenvalues=eigen(C,symmetric=TRUE)$values)
}

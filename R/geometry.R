# Floating-point evaluations of the three examples; exact signs are separate.
run_geometry <- function() {
  d<-2.8; r<-.68; beta<-.99; tau<-.01
  means<-rbind(c(0,0),c(d,0),c(d/2,sqrt(3)*d/2))
  A<-exp(-d*d/2+r*d); c0<-c(.51,0)
  L<-matrix(c(2.72,2.7,2.7,2.7),2)
  dc<-c(sqrt(3)*d/2,d/2)-c0
  triangle<-list(means=means,weights=rep(1/3,3),covariance=diag(2),radius=r,
    odds_cap=A,first_matrix_slack=min(beta-c0[1]^2-2*A*2.72,beta-2*A*2.7),
    multiplier_determinant=det(L),second_matrix_determinant=det(beta*diag(2)+L-tcrossprod(dc)),
    separation_slack=d*d-(3+tau+r*r+2/r^2),
    saturated=saturated_moments(means,rep(1/3,3),1,r),
    minimum_radius=sqrt(rho_squared(d,2,tau)),
    trace_at_minimum_radius=saturated_moments(means,rep(1/3,3),1,sqrt(rho_squared(d,2,tau)))$trace,
    analytic_sdp_optimum=2*A*(sqrt(3)*d/2)^2/(1+2*A)^2,
    mode_search=locate_modes(means))
  v<-matrix(c(0,3,7),ncol=1); w<-c(50,50,1)/101
  radii<-sqrt(rho_squared(c(3,3,4),1,.01))
  local<-lapply(1:3,function(j)saturated_moments(v,w,j,radii[j]))
  common<-lapply(1:3,function(j)saturated_moments(v,w,j,max(radii)))
  # A conservative derivative bound verifies that the saturated odds maximize
  # variance for this example. This property is not assumed for general models.
  derivative_bounds<-lapply(1:3,function(j) {
    delta<-as.numeric(v[-j]-v[j]); a<-common[[j]]$odds
    (abs(delta)-sum(a*abs(delta)))^2-sum(a*delta^2)
  })
  stopifnot(min(unlist(derivative_bounds))>0)
  radius_table<-data.frame(component=1:3,local_radius=radii,
    local_variance=vapply(local,`[[`,0,"trace"),common_variance=vapply(common,`[[`,0,"trace"))
  v2<-rbind(c(0,0),c(2.3,0),c(0,20)); rs<-c(.6,.6,.06); nu<-1.001
  Ms<-list(diag(c(1,0)),diag(c(1,0)),diag(c(0,1)))
  eps<-vapply(1:3,function(j)sum(vapply(setdiff(1:3,j),function(h) {
    delta<-v2[h,]-v2[j,]; dd<-sqrt(sum(delta^2))
    exp(-dd^2/2+sqrt(2)*dd)*drop(t(delta)%*%(diag(2)-Ms[[j]])%*%delta)
  },0)),0)
  budget<-list(means=v2,radii=rs,nu=rep(nu,3),beta=.95,tau=.01,residual_sums=eps,
    trace_budget_slacks=nu-1-eps/.95,
    separation_slacks=c(2.3,2.3,20)^2-(nu+1+.01+rs^2+nu/rs^2),
    saturated_traces=vapply(1:3,function(j)saturated_moments(v2,rep(1/3,3),j,rs[j])$trace,0),
    full_rank_floor_slack=2.3^2-6.01)
  stopifnot(triangle$first_matrix_slack>0,triangle$separation_slack>0,
    triangle$trace_at_minimum_radius>1, all(radius_table$local_variance<.9),
    radius_table$common_variance[3]>1,all(budget$separation_slacks>0),all(eps>0))
  list(triangle=triangle,local_radius=list(table=radius_table,variance_derivative_lower_bounds=derivative_bounds),
       local_budget=budget,arithmetic="IEEE double evaluations; exact sign certificates are in exact_geometry.json")
}

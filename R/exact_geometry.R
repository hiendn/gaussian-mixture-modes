# Exact rational interval checks for the three geometric manuscript examples.
# All decisions below use gmp rationals. Decimal conversions are for display.
# The original Python calculations are preserved in original/code/.

ex_q <- function(x) {
  if (inherits(x, "bigq")) return(x)
  if (inherits(x, "bigz")) return(gmp::as.bigq(x))
  if (length(x) != 1L) stop("Expected one exact scalar")
  if (is.numeric(x) && (!is.finite(x) || x != trunc(x)))
    stop("Use a quoted decimal or fraction for noninteger exact inputs")
  x <- as.character(x)
  if (grepl("^[+-]?[0-9]+\\.[0-9]+$", x)) {
    n <- nchar(sub(".*\\.", "", x))
    digits <- sub("^[+-]", "", gsub(".", "", x, fixed = TRUE))
    digits <- sub("^0+", "", digits)
    if (!nzchar(digits)) digits <- "0"
    if (startsWith(x, "-")) digits <- paste0("-", digits)
    return(gmp::as.bigq(digits) / gmp::as.bigz(10)^n)
  }
  gmp::as.bigq(x)
}

ex_i <- function(lo, hi = lo) {
  if (inherits(lo, "exact_interval") && missing(hi)) return(lo)
  lo <- ex_q(lo); hi <- ex_q(hi)
  if (!isTRUE(lo <= hi)) stop("Reversed interval")
  structure(list(lo = lo, hi = hi), class = "exact_interval")
}

Ops.exact_interval <- function(e1, e2) {
  a <- ex_i(e1)
  if (missing(e2)) {
    if (.Generic == "-") return(ex_i(-a$hi, -a$lo))
    if (.Generic == "+") return(a)
    stop("Unsupported unary interval operation")
  }
  if (.Generic == "^") {
    if (length(e2) != 1L || e2 != 2L) stop("Only interval squares are required")
    hi <- max(a$lo^2, a$hi^2)
    lo <- if (a$lo <= 0 && a$hi >= 0) ex_q(0) else min(a$lo^2, a$hi^2)
    return(ex_i(lo, hi))
  }
  b <- ex_i(e2)
  switch(.Generic,
    "+" = ex_i(a$lo + b$lo, a$hi + b$hi),
    "-" = ex_i(a$lo - b$hi, a$hi - b$lo),
    "*" = {
      endpoints <- c(a$lo*b$lo, a$lo*b$hi, a$hi*b$lo, a$hi*b$hi)
      ex_i(min(endpoints), max(endpoints))
    },
    "/" = {
      if (b$lo <= 0 && b$hi >= 0) stop("Interval division through zero")
      a * ex_i(1/b$hi, 1/b$lo)
    }, stop("Unsupported interval operation: ", .Generic))
}

ex_sum <- function(x) Reduce(`+`, x, init = ex_i(0))

ex_integer_sqrt <- function(n) {
  n <- gmp::as.bigz(n)
  if (n < 0) stop("Negative square-root argument")
  if (n < 2) return(n)
  x <- n
  repeat {
    y <- (x + n %/% x) %/% 2
    if (y >= x) return(x)
    x <- y
  }
}

ex_sqrt_q <- function(x, digits = 50L) {
  x <- ex_q(x)
  if (x < 0) stop("Negative square-root argument")
  den <- gmp::denominator(x); num <- gmp::numerator(x)
  scale <- gmp::as.bigz(10)^digits
  target <- num*den*scale^2
  k <- ex_integer_sqrt(target)
  if (!(k^2 <= target && (k+1)^2 > target)) stop("Integer root enclosure failed")
  upper <- if (k^2 == target) k else k+1
  ex_i(gmp::as.bigq(k)/(den*scale), gmp::as.bigq(upper)/(den*scale))
}

ex_sqrt <- function(x, digits = 50L) {
  x <- ex_i(x)
  ex_i(ex_sqrt_q(x$lo, digits)$lo, ex_sqrt_q(x$hi, digits)$hi)
}

ex_exp_negative_q <- function(x, degree = 100L) {
  x <- ex_q(x)
  if (!(x >= 0 && x < degree+2)) stop("Taylor-tail precondition failed")
  term <- total <- ex_q(1)
  for (k in seq_len(degree)) {
    term <- term*x/k
    total <- total+term
  }
  next_term <- term*x/(degree+1)
  upper <- total+next_term/(1-x/(degree+2))
  ex_i(1/upper, 1/total)
}

ex_exp_negative <- function(x, degree = 100L) {
  x <- ex_i(x)
  # exp(-x) is decreasing. Each endpoint is enclosed independently.
  ex_i(ex_exp_negative_q(x$hi, degree)$lo,
       ex_exp_negative_q(x$lo, degree)$hi)
}

ex_bound <- function(name, quantity, direction, threshold = 0) {
  x <- ex_i(quantity); threshold <- ex_q(threshold)
  passed <- switch(direction,
    less_than = x$hi < threshold,
    greater_than = x$lo > threshold,
    stop("Unknown strict relation"))
  if (!isTRUE(passed)) stop("Exact check failed: ", name)
  list(quantity = name, certified_relation = direction,
       rational_threshold = as.character(threshold),
       rational_lower_endpoint = as.character(x$lo),
       rational_upper_endpoint = as.character(x$hi),
       threshold_decimal_for_display = as.numeric(threshold),
       midpoint_for_display_only = as.numeric((x$lo+x$hi)/2),
       sign_decision_arithmetic = "exact rational", passed = TRUE)
}

ex_separation <- function(d2, r, nu, tau) {
  ex_i(d2)-(ex_i(nu)+1+ex_i(tau)+ex_i(r^2)+ex_i(nu/r^2))
}

exact_local_radius <- function() {
  q <- ex_q
  means <- lapply(c(0,3,7), q)
  weights <- list(q("50/101"),q("50/101"),q("1/101"))
  radii <- lapply(c("0.383","0.383","0.269"), q)
  distances <- lapply(c(3,3,4), q)
  checks <- list()
  for (j in seq_len(3L)) {
    others <- setdiff(seq_len(3L),j)
    deltas <- lapply(others, function(h) means[[h]]-means[[j]])
    odds <- lapply(seq_along(others), function(k) {
      h <- others[k]; d <- abs(deltas[[k]])
      ex_i(weights[[h]]/weights[[j]]) * ex_exp_negative(d^2/2-radii[[j]]*d)
    })
    S <- ex_sum(odds); denominator <- 1+S
    centre <- ex_sum(lapply(seq_along(odds), function(k) odds[[k]]*ex_i(deltas[[k]])))/denominator
    second <- ex_sum(lapply(seq_along(odds), function(k) odds[[k]]*ex_i(deltas[[k]]^2)))/denominator
    variance <- second-centre^2
    checks <- c(checks, list(
      ex_bound(paste0("j",j,"_separation_slack"),
        ex_separation(distances[[j]]^2,radii[[j]],q(1),q("0.01")),
        "greater_than",c("0.026","0.026","0.098")[j]),
      ex_bound(paste0("j",j,"_S"),S,"less_than","0.05"),
      ex_bound(paste0("j",j,"_saturated_trace"),variance,"less_than",c("0.30","0.30","0.72")[j])))
  }
  lower <- q("0.382")
  checks <- c(checks,list(ex_bound("common_radius_.382_separation_failure",
    -ex_separation(q(9),lower,q(1),q("0.01")),"greater_than","0.0088")))
  a <- 50*ex_exp_negative(q(8)-4*lower)
  checks <- c(checks,list(ex_bound("fixed_two_point_robust_variance_in_common_ball",
    16*a/(1+a)^2,"greater_than","1.06")))
  list(means=c(0,3,7),weights=c("50/101","50/101","1/101"),
       beta="0.9",tau="0.01",radii=c("0.383","0.383","0.269"),
       exponential_enclosure="Rational Taylor sum through degree 100 with geometric positive-tail bound",
       checks=checks)
}

exact_local_budget <- function() {
  q <- ex_q
  beta <- q("0.95"); tau <- q("0.01"); nu <- q("1.001")
  radii <- lapply(c("0.6","0.6","0.06"),q)
  d2 <- list(q("5.29"),q("5.29"),q(400))
  eps <- q(400)/gmp::as.bigz(2)^170
  checks <- list(
    ex_bound("coarse_residual_epsilon_upper",eps,"less_than",q(1)/gmp::as.bigz(10)^45),
    ex_bound("effective_budget_slack_lower",nu-1-eps/beta,"greater_than","0.0009"))
  for (j in seq_len(3L)) checks <- c(checks,list(ex_bound(paste0("j",j,"_separation_slack"),
    ex_separation(d2[[j]],radii[[j]],nu,tau),"greater_than",c("0.138","0.138","119.9")[j])))
  x <- q("1.265"); term <- total <- q(1)
  for (k in seq_len(8L)) {term <- term*x/k; total <- total+term}
  stopifnot(total > q("1000/283"))
  pair_odds <- q("0.283")
  pair_var <- q("5.29")*pair_odds/(1+pair_odds)^2
  remote_odds <- q(1)/gmp::as.bigz(2)^188
  isolated_odds <- q(1)/gmp::as.bigz(2)^198
  checks <- c(checks,list(
    ex_bound("close_pair_variance_upper",pair_var,"less_than","0.91"),
    ex_bound("j1_j2_saturated_trace_upper",pair_var+576*remote_odds,"less_than","0.911"),
    ex_bound("j1_j2_S_upper",pair_odds+remote_odds,"less_than","0.284"),
    ex_bound("j3_saturated_trace_upper",2*441*isolated_odds,"less_than",q(1)/gmp::as.bigz(10)^50),
    ex_bound("j3_S_upper",2*isolated_odds,"less_than",q(1)/gmp::as.bigz(10)^50),
    ex_bound("old_global_q2_floor_failure",q("6.01")-q("5.29"),"greater_than","0.7")))
  list(means=list(c(0,0),c("2.3","0"),c(0,20)),weights=rep("1/3",3),
       beta="0.95",tau="0.01",nu="1.001",radii=c("0.6","0.6","0.06"),
       tail_proof="sqrt(2)<3/2, monotonicity of -d^2/2+s*d for d>=20, and e>2",
       checks=checks)
}

exact_triangle <- function() {
  q <- ex_q
  d <- q("2.8"); r <- q("0.68"); beta <- q("0.99"); tau <- q("0.01"); nu <- q(2)
  A <- ex_exp_negative(d^2/2+r*(-d))
  a <- ex_sqrt(q(3))*ex_i(d)/2; b <- d/2; c <- q("0.51")
  u <- q("2.72"); v <- q("2.7"); z <- q("2.7")
  H00 <- ex_i(beta+u)-(a-ex_i(c))^2; H11 <- beta+v-b^2; H01 <- ex_i(z)-(a-ex_i(c))*ex_i(b)
  quantities <- list(
    separation_slack=ex_separation(d^2,r,nu,tau),
    invariance_slack=ex_i(r)-2*A*ex_i(d-r),
    lambda_determinant=ex_i(u*v-z^2),
    main_slack_x=ex_i(beta-c^2)-2*A*ex_i(u),
    main_slack_y=ex_i(beta)-2*A*ex_i(v),
    second_slack_x=H00,
    second_slack_y=ex_i(H11),
    second_determinant=H00*ex_i(H11)-H01^2)
  thresholds <- c("0.04234","0","0.05399","0.00536","0","0","0","0.07448")
  checks <- lapply(seq_along(quantities),function(k)
    ex_bound(names(quantities)[k],quantities[[k]],"greater_than",thresholds[k]))
  root_a <- d^2-nu-1-tau
  rmin <- ex_sqrt(2*ex_i(nu)/(ex_i(root_a)+ex_sqrt(root_a^2-4*nu)))
  Amin <- ex_exp_negative(ex_i(d^2/2)-ex_i(d)*rmin)
  trace <- 2*Amin*ex_i(d^2)/(1+2*Amin)-3*Amin^2*ex_i(d^2)/(1+2*Amin)^2
  checks <- c(checks,list(ex_bound("trace_at_minimum_admissible_radius",trace,"greater_than",1)))
  list(side_length="2.8",radius="0.68",beta="0.99",tau="0.01",nu="2",
       square_root_enclosure="Exact integer square root after rational scaling by 10^50",
       exponential_enclosure="Rational Taylor sum through degree 100 with geometric positive-tail bound",
       checks=checks,
       minimum_radius_decimal_for_display=as.numeric((rmin$lo+rmin$hi)/2),
       trace_at_minimum_radius_decimal_for_display=as.numeric((trace$lo+trace$hi)/2))
}

run_exact_geometry <- function() {
  if (!requireNamespace("gmp",quietly=TRUE)) stop("Install the gmp package for exact checks")
  triangle <- exact_triangle(); radii <- exact_local_radius(); budgets <- exact_local_budget()
  n <- length(triangle$checks)+length(radii$checks)+length(budgets$checks)
  list(arithmetic="Native R with gmp exact rationals throughout all strict sign decisions",
       scope="Three manuscript geometric examples; floating-point fitting iterates are not interval certified",
       gmp_version=as.character(utils::packageVersion("gmp")),
       checks_passed=n,triangle_matrix_example=triangle,
       one_dimensional_local_radius_example=radii,
       two_dimensional_local_budget_example=budgets)
}

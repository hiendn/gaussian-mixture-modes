#!/usr/bin/env Rscript
# Run from the repository root. Installs into this project only.
if(getRversion()<"4.4.0") stop("Use R 4.4 or later for the current CVXR dependencies.")
lib <- file.path(getwd(),".R-library")
dir.create(lib,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"))
needed <- c("jsonlite","gmp","CVXR","clarabel")
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
for(p in intersect(c("CVXR","clarabel"),setdiff(needed,missing))) {
  minimum<-if(p=="CVXR")"1.8.1" else "0.11.1"
  if(packageVersion(p)<minimum) missing<-union(missing,p)
}
if(length(missing)) install.packages(missing,lib=lib)
cat("R library:",lib,"\n")
for(p in needed) cat(p,as.character(packageVersion(p)),"\n")
